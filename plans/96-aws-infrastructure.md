# 96 — AWS for storage and burst compute: a Terraform-managed, cost-gated cloud leg

Status: 📋 planned 2026-09-26 — **no infrastructure exists yet and nothing here has been
applied.** Phase 0 (account guardrails) is **done and lives in the
[`aws-cloud`](https://github.com/joebarbere/aws-cloud) repo** as of 2026-09-26; phases 1–3
(this repo's own bucket and compute) have not started. Prices are us-east-1 list prices pulled from AWS's public price-list files on
2026-09-26 (S3 offer `publicationDate` 2026-09-26T01:55Z; AWSDataTransfer offer); **spot prices
are third-party averages** (instances.vantage.sh) and move daily — re-check before any run.

## Context

Two local constraints are real, and one is not:

1. **Disk (real).** 70 GB free of 929 GB on 2026-09-26 (93% full). `ideas.md` was written
   assuming ~275 GB. Apertif TD DR2 (~0.8 PB), LoTSS DR3 images, and any cutout-heavy sweep do
   not fit.
2. **CUDA-only tooling (real).** The local RX 7600 XT runs pure PyTorch on ROCm; Heimdall,
   AstroAccelerate, dedisp, CuPy- and RAPIDS-locked tools do not run. Those are exactly the
   independent oracles `torchfdmt`/`singlepulse` lack (ideas.md L6), and the `torchfdmt`
   paper's "device-portable" claim has only ever run on CPU + ROCm.
3. **Raw throughput (not real, mostly).** Most open plans are CASDA-bound (the network is the
   bottleneck, not the CPU) or fit the 16 GB local GPU. Renting compute to go faster is rarely
   worth it here; renting it to go *somewhere the local machine cannot* is.

There is also a second motive: Joe is studying for **AWS Certified AI Practitioner (AIF-C01)**
in the sibling `aws-ai` repo (`~/dev/github/joebarbere/aws-ai`), which already has Identity
Center SSO (`joebarbere-admin` profile), a cost-seatbelt deny policy, and Terraform modules with
`enable_*` cost gates. This plan reuses those conventions rather than inventing new ones. Exam
practice stays in `aws-ai`; this repo only provisions what the research needs.

## Recommendation (the opinion asked for)

**Yes to a Terraform directory, kept small, and yes to AWS — for storage and CUDA, not for
speed.** Specifically:

- `infra/terraform/` (not top-level `terraform/`): it is not part of the Python package, and
  `infra/` leaves room for a Packer/AMI or container definition later.
- **One shared account, separated by tags (owner decision 2026-09-26).** `default_tags`
  sets `Project = jansky-research` here and `Project = aws-ai` in the study repo; the two
  projects keep **separate Terraform state**, so a `terraform destroy` in one can only remove
  what that state created. What tags do *not* give, and how the plan covers it:
  - *Cost split* — only after `Project` is **activated** as a cost-allocation tag, and only
    from activation onward (not retroactive). Untagged spend (e.g. some data-transfer lines,
    anything created by hand) lands in "no tag", so the account-wide budget is the backstop.
  - *Blast radius* — tags don't limit permissions. The seatbelt is account-wide, so it
    constrains both projects; that is intended (an allowlist of cheap instance types suits the
    study repo too). Name prefixes (`jansky-`) keep resources distinguishable in the console.
  - A separate member account remains the upgrade path if spend grows or data becomes
    something that must never be deleted by accident.
- **The science never depends on AWS.** Same rule as the GPU: every cloud result must reproduce
  at small scale on the CPU path, and no committed number may depend on which device or
  provider produced it.

## Prices that drive the design (us-east-1, 2026-09-26)

**Storage (S3, per GB-month; list prices from the offer file)**

| Class | $/GB-mo | ≈ $/TB-mo | Use |
|---|---|---|---|
| Standard (first 50 TB) | 0.023 | 23.55 | working set during a run |
| Standard-IA | 0.0125 | 12.80 | downloaded archive products kept for re-runs |
| Glacier Instant Retrieval | 0.004 | 4.10 | cold raw data you might need again |
| Glacier Flexible Retrieval | 0.0036 | 3.69 | archival only |

IA and Glacier classes add per-GB retrieval fees and minimum storage durations — check both
before choosing a lifecycle rule. EBS gp3 is roughly $0.08/GB-mo (standard list price, **not**
re-pulled here), i.e. ~3.5× S3 Standard: keep bulk data in S3, not on volumes.

**Egress is the trap.** Data transfer out to the internet: **$0.09/GB** for the first 10 TB/month
after the 100 GB/month free tier. Bringing 1 TB home costs ~$83. Ingress is free. So: **download
archives straight into AWS, compute there, bring back only results** (JSON/CSV/figures are MB).

**Compute (EC2 Linux, $/hour)**

| Instance | vCPU / RAM | GPU | On-demand | Spot (avg) | Note |
|---|---|---|---|---|---|
| g4dn.xlarge | 4 / 16 GB | T4 16 GB | 0.526 | ~0.26 | cheapest CUDA; fine for oracle cross-checks |
| g5.xlarge | 4 / 16 GB | A10G 24 GB | 1.006 | ~0.50 | best spot discount (~50%) |
| g6.xlarge | 4 / 16 GB | L4 24 GB | 0.805 | ~0.62 | cheapest on-demand 24 GB |
| g6e.xlarge | 4 / 32 GB | L40S 48 GB | 1.861 | ~1.82 | only if a model needs >24 GB |
| c7i.4xlarge | 16 / 32 GB | — | 0.714 | ~0.26 | CPU sweeps |
| c7a.8xlarge | 32 / 64 GB | — | 1.642 | ~0.61 | CPU sweeps, more cores |
| p4d / p5 | 8× A100 / 8× H100 | — | 21.96 / 55.04 | — | **out of scope**; the seatbelt denies them |

The local RX 7600 XT is 16 GB; every 24 GB card above is a step up in memory, not only in CUDA.

## Worked costs for the concrete uses

| Use | Assumptions | Estimate |
|---|---|---|
| **Idle stack** | budget + empty bucket + roles; no instances | ~$0–1/month |
| **Storage offload** | 1 TB of archive products in Standard-IA | ~$13/month |
| **CUDA validation of `torchfdmt`** (phase 2) | g6.xlarge on-demand, ~3 h incl. setup, repo's Crab file | ~$2.50 |
| **Plan 51 open-weights classifier, CUDA baseline** | g5.xlarge spot, ~50 h | ~$25 |
| **Plan 61 Apertif subset + Heimdall oracle (L6)** | ~2 TB staged in S3 Standard for a month, g6.xlarge spot ~200 h | ~$47 storage + ~$125 compute ≈ **$170 one-off** |
| **Mistake: bring that 2 TB home** | egress | **~$175** — as much as the run |
| **Mistake: forget a g6.xlarge on-demand** | 30 days | **~$580** |
| **Mistake: a NAT gateway** | hourly + per-GB processing | tens of $/month, forever |

The last three rows are why the guardrails in phase 0 come before anything else.

## Terraform layout (mirrors `aws-ai`)

**Three repos, three states, one account.** Account-wide resources — the seatbelt, the
cost-allocation tags, every budget (account and per-project) and the anomaly monitor — are owned
by [`aws-cloud`](https://github.com/joebarbere/aws-cloud) (decided 2026-09-26: two states
managing one budget or one permission-set policy would each undo the other's apply). This repo's
Terraform holds only jansky-research's own resources:

```
infra/terraform/
  versions.tf  providers.tf   # aws ~> 5.70; default_tags Project=jansky-research, ManagedBy=terraform
  variables.tf                # region, name_prefix, enable_* gates (budgets live in aws-cloud)
  main.tf  outputs.tf
  terraform.tfvars.example    # tfvars + state gitignored, exactly as in aws-ai
  modules/
    storage/     always on    one private bucket: public-access block, SSE-S3, lifecycle
                              (raw/ → IA at 30 d → Glacier IR at 90 d; results/ stays Standard;
                              abort incomplete multipart uploads at 7 d), S3 gateway VPC endpoint (free)
    compute/     gated off    launch template for spot GPU/CPU: instance profile scoped to the
                              bucket, SSM Session Manager (no SSH key, no inbound rules), default
                              VPC public subnet (NO NAT gateway), gp3 root, IMDSv2 required,
                              shutdown_behavior = terminate + a hard max-lifetime in user-data
    sagemaker/   gated off    execution role + bucket access for managed-spot training jobs only
                              (plan 51); no domain, no endpoints — those bill hourly
```

Any instance type a job needs must be on the seatbelt's EC2 allowlist in `aws-cloud` first.

Local state (gitignored) is fine for one operator; an S3 backend with a lock is a later
upgrade, not a prerequisite. CI gets `terraform fmt -check` + `terraform validate` (no
credentials needed); **never** `plan`/`apply` in CI.

How jobs run: user-data installs `uv`, clones this repo **at a tag**, runs one named slice
command with `--out` pointed at a scratch dir, syncs results to `s3://…/results/<slice>/<run>/`,
then terminates. Spot interruptions are absorbed by the resumable per-row CSV pattern the
real-leg scripts already use (`load_done` in `scripts/stokesv_discovery_real.py`). Results come
home with `aws s3 sync` and are committed as evidence exactly like local runs — with the
`benchmark_device` / `benchmark_hardware` fields recording where they ran.

## Phases

0. **Guardrails in the shared account — $0.** State on 2026-09-26 (read-only audit): the
   `Project`, `ManagedBy` and `Purpose` tags exist but are **Inactive** for cost allocation;
   **no budgets** and **no cost-anomaly monitors** exist; the `AdministratorAccess` permission
   set has **no inline policy**, i.e. `aws-ai`'s `cost-seatbelt.json` is written but not
   attached. No EC2 instances were running. Steps: activate `Project` (and `ManagedBy`) as
   cost-allocation tags; an account-wide monthly budget plus one `Project`-filtered budget per
   repo; a cost-anomaly monitor; attach a merged seatbelt (aws-ai's statements + the EC2
   allowlist and NAT/p4d/p5 denies) to the permission set. Test the seatbelt by trying to
   launch a denied instance type and confirming the deny.
   **Done 2026-09-26 (steps 1–3, via the AWS CLI, owner-approved):** `Project` and `ManagedBy`
   activated as cost-allocation tags (split visible from this date only); budgets
   `account-monthly-25` ($25/month, account-wide) plus `jansky-research-monthly` and
   `aws-ai-monthly` ($25 each, filtered on `user:Project`) — **alerts only**, email at 50/80/100%
   actual and 100% forecast; cost-anomaly monitor `services-anomaly-monitor` (per service) with a
   daily email for anomalies ≥ $5.
   **Seatbelt attached 2026-09-26 (owner-approved):** the merged seatbelt became the inline policy
   on the `AdministratorAccess` permission set, provisioned to the account. It is `aws-ai`'s
   `iam/cost-seatbelt.json` (region lock to us-east-1; no Kendra / OpenSearch Serverless /
   SageMaker endpoints / RDS) **plus** an EC2 instance-type allowlist and denies on NAT gateways,
   fleets, legacy spot requests, dedicated hosts, reservations and Savings Plans. **The attached
   policy is the source of truth; `aws-ai`'s file is now a subset of it.** Access Analyzer: no
   findings. Proven by dry-run, not assumed: `p5.48xlarge` and `g6e.2xlarge` → denied;
   `g6.xlarge`, `g5.xlarge` spot, `t3.micro` → allowed; `create-nat-gateway` on a real
   default-VPC subnet → denied; `request-spot-instances` → denied; any call in `us-west-2` →
   denied; STS and Budgets still work. Spot runs go through `run-instances
   --instance-market-options`, which the allowlist governs.
   **Imported into Terraform 2026-09-26 — the [`aws-cloud`](https://github.com/joebarbere/aws-cloud)
   repo now owns all of the above** (8 resources; the post-import apply added only default tags
   and relabelled the policy's `Id`; a second plan showed no changes; denies re-proven by
   dry-run). The CLI copies, `infra/seatbelt.json`, `aws-ai/iam/cost-seatbelt.json` and the
   `make seatbelt-check` script were removed: `terraform plan` in `aws-cloud` is the drift check.
   To change the seatbelt or a budget, change it there.
1. **Storage — ~$0 until used.** `storage` module; upload nothing yet. Verify the public-access
   block and lifecycle rules with `aws s3api get-bucket-*`.
2. **One CUDA validation run — ~$2.50.** `compute` module, g6.xlarge on-demand, the `torchfdmt`
   Crab recover-a-known with `--device cuda`. Success = the same DM to the committed tolerance
   and a CPU-oracle parity pass on NVIDIA. This strengthens an existing paper claim
   ("device-portable") at trivial cost, and proves the whole pipeline end to end.
   `terraform destroy` afterwards; `terraform state list` must be empty of compute.
   **Done 2026-09-26 — PASS, $0.29.** Evidence: `results/cuda_validation_2026-09-26.json`.
   - Ran on **g5.xlarge (NVIDIA A10G 24 GB)**, not g6.xlarge: the L4 had no on-demand capacity
     in any us-east-1 zone at launch (the seatbelt's region lock forbids trying elsewhere, as
     intended). Launched from the `jansky-gpu` template (`infra/terraform/`), job
     `infra/jobs/cuda_validation.sh` via SSM Run Command, results fetched over SSM (no bucket
     needed for KB-scale output), terminated by hand after the copy-home. 17 min 11 s of
     instance time = **$0.288** at $1.006/h, + ~$0.002 of gp3; the plan estimated $2.50.
   - torch **2.13.0+cu130** (the `uv.lock` version; the committed ROCm run was 2.12.1).
     2.13.0 has no `cu128` wheel — the first attempt failed on that, a second on a relative
     output path in my job script; both were setup failures before any science ran.
   - **36/36** fdmt/singlepulse/torchdsp tests pass. FDMT plane CPU vs CUDA: **bit-identical**
     (max diff 0.0). Brute dedispersion: max relative diff 1.8e-7 (float summation order),
     same best DM.
   - **The Crab recover-a-known with the science leg on CUDA reproduces every committed CPU
     number exactly** (real DM 56.59, S/N 14.0, position 3501, 2.9 trials, p_pos 0.00353, null
     5.23/6.03/6.26, p_null 0.0299; synthetic 56.63 / 33.379 ms).
   - Bonus: `torchdsp`'s cross-device kernel check on NVIDIA is **identical field-for-field**
     to the committed AMD one, so that paper's portability claim now holds across two GPU
     vendors.
   - A10G timings (not comparable to the committed AMD rows; different host): brute 0.79 s
     vs 85.4 s on the 4-vCPU host; FDMT 0.84 s GPU vs 0.68 s CPU (FDMT is already cheap);
     SumThreshold is *slower* on the GPU (12.0 s vs 4.8 s) — its iteration structure does not
     parallelise well, consistent with the torchdsp paper's framing.
   - Terraform kept: the role/SG/launch template cost nothing idle; `terraform destroy` in
     `infra/terraform/` removes them.
3. **First real workload** — plan 61 (Apertif subset + Heimdall oracle) or plan 51 (classifier
   CUDA baseline), whichever is picked first. Budget the run in advance from the table above.

## Verification

- `terraform validate` + `fmt -check` clean; `tflint` optional.
- Seatbelt proven by a denied call, not assumed (the CLAUDE.md lesson: test a guard through the
  path that actually runs it).
- Phase 2 reproduces the committed `torchfdmt` Crab DM; results JSON records
  `benchmark_hardware` = the NVIDIA card, and the science numbers are unchanged.
- After every session: `terraform state list` shows no compute; Cost Explorer shows the run's
  spend within the estimate. Record actual vs estimated cost in this file.

## Owner decisions needed before phase 0

- ~~Separate member account or shared?~~ **Decided 2026-09-26: shared account, tag-separated.**
- ~~Budget and alert mode~~ **Decided 2026-09-26: $25/month account-wide, alerts only.**
- ~~Attach the merged seatbelt?~~ **Done 2026-09-26** (see phase 0).
- ~~Where do account-wide resources live?~~ **Decided 2026-09-26: a separate `aws-cloud` repo.**
