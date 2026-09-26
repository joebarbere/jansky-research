# infra/ — this repo's AWS resources (plan 96)

The plan, prices and phases are in [`plans/96-aws-infrastructure.md`](../plans/96-aws-infrastructure.md).

## Three repos, one account

| Repo | Owns |
|---|---|
| [`aws-cloud`](https://github.com/joebarbere/aws-cloud) | **everything account-wide**: the cost seatbelt (explicit-deny inline policy on `AdministratorAccess`), cost-allocation tag activation, all budgets (account-wide $25/month and per-`Project`), the cost-anomaly monitor |
| this repo, `infra/terraform/` (not built yet) | jansky-research's own bucket and spot compute, tagged `Project = jansky-research` |
| [`aws-ai`](https://github.com/joebarbere/aws-ai) | the AIF-C01 study modules |

Each keeps its own Terraform state, so a `destroy` in one can only remove what that state
created. **Do not declare budgets, tag activations or IAM guardrails here** — two states
managing one account-wide object undo each other's applies.

## What that means when working here

- **Instance types:** the seatbelt denies `ec2:RunInstances` for any type not on its allowlist
  (small, CPU, and single-GPU `g4dn`/`g5`/`g6` types). A job needing another type means a change
  to `aws-cloud/terraform/seatbelt.json` first, then `make plan && make apply` there, then a
  `--dry-run` to prove it.
- **Spot:** use `run-instances --instance-market-options MarketType=spot`; the fleet and legacy
  spot-request APIs are denied. No NAT gateways; only `us-east-1`.
- **Cost:** tag everything `Project = jansky-research` (provider `default_tags`), or the
  `jansky-research-monthly` budget cannot see it. The account-wide budget catches the rest.
- **Drift:** `make drift` in `aws-cloud` (exit 0 = the account matches the code).

## What exists (2026-09-26)

`infra/terraform/` — free until an instance runs: an EC2 role with SSM access only, a security
group with **no inbound rules**, and the `jansky-gpu` launch template (NVIDIA-driver base AMI,
IMDSv2, encrypted gp3 root, shutdown = terminate, user-data schedules a shutdown
`max_lifetime_minutes` after boot). State is local and gitignored.

One-shot GPU job, as run for plan 96 phase 2:

```fish
set -x AWS_PROFILE joebarbere-admin
terraform -chdir=infra/terraform init; and terraform -chdir=infra/terraform apply
aws ec2 run-instances --launch-template LaunchTemplateName=jansky-gpu \
    --instance-type g5.xlarge --subnet-id <default subnet in an AZ with capacity>
# wait for SSM PingStatus Online, then:
aws ssm send-command --instance-ids <id> --document-name AWS-RunShellScript \
    --parameters 'commands=["curl -fsSL https://raw.githubusercontent.com/joebarbere/jansky-research/<sha>/infra/jobs/cuda_validation.sh -o /opt/cv.sh && bash /opt/cv.sh <sha>"]'
# fetch /opt/job/out/* over SSM, then terminate the instance yourself
```

- **Capacity is not guaranteed.** g6.xlarge had no on-demand capacity in any us-east-1 zone on
  2026-09-26; loop over subnets and fall back to another allowlisted type.
- **The job pins a commit**, so push before launching; results go to a scratch dir on the
  instance, never into a `results/` tree.
- **Terminate after copying results home** rather than waiting for the max-lifetime shutdown.

`infra/jobs/cuda_validation.sh` — the torch-fdmt Crab recover-a-known + CPU/CUDA parity +
torch-dsp cross-device checks; its evidence is `results/cuda_validation_2026-09-26.json`.
