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

Nothing is provisioned for this repo yet.
