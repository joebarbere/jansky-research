# infra/ — AWS for the cloud leg (plan 96)

What lives here, and what is live in the account. The plan, prices and phases are in
[`plans/96-aws-infrastructure.md`](../plans/96-aws-infrastructure.md).

## `seatbelt.json` — the source of truth for the account's cost seatbelt

This repo and the sibling [`aws-ai`](https://github.com/joebarbere/aws-ai) study repo share **one
AWS account**, separated by the `Project` cost-allocation tag and by separate Terraform state. A
permission set holds one inline policy, so there is **one** seatbelt for both, and this file is it.
It has been the inline policy on the `AdministratorAccess` permission set since 2026-09-26.

`aws-ai/iam/cost-seatbelt.json` is a **byte-for-byte copy**. Edit here, never there.

| Statement | Denies |
|---|---|
| `RegionLock` | everything outside `us-east-1` (global services exempt) |
| `NoHourlyBilledResourcesWithoutEditingThisPolicy` | Kendra indexes, OpenSearch Serverless, SageMaker endpoints, RDS/Aurora |
| `Ec2InstanceTypeAllowlist` | `ec2:RunInstances` for any type outside: `t3.micro`/`t3.small`/`t4g.small`, `c7i.4xlarge`, `c7a.4xlarge`/`8xlarge`, `m7a.2xlarge`, `r7a.2xlarge`, `g4dn.xlarge`, `g5.xlarge`/`2xlarge`, `g6.xlarge`/`2xlarge`, `g6e.xlarge` |
| `NoNatGatewayNoFleetsNoCommitments` | NAT gateways, `CreateFleet`, spot fleets, legacy spot requests, dedicated hosts, capacity reservations, Reserved Instances, Savings Plans |

Spot instances go through `run-instances --instance-market-options MarketType=spot`, which the
allowlist governs; the fleet APIs are denied because they launch through paths it does not see.

**Comments:** IAM policy JSON allows none, and IAM rejects unknown keys. The optional top-level
`Id` is the one free-text field, so it names the source of truth and the copy. IAM ignores `Id`
when evaluating (Access Analyzer: no findings; IAM accepted it on attach).

### Changing it

1. Edit `infra/seatbelt.json`; `aws accessanalyzer validate-policy --policy-type IDENTITY_POLICY
   --policy-document file://infra/seatbelt.json`.
2. Copy it to `../aws-ai/iam/cost-seatbelt.json` and commit in both repos.
3. `aws sso-admin put-inline-policy-to-permission-set …` then `provision-permission-set …
   --target-type ALL_PROVISIONED_ACCOUNTS` (an edited permission set does nothing until
   re-provisioned).
4. `make seatbelt-check` — compares the aws-ai copy **and the live attached policy** with this
   file, and exits non-zero on drift. Needs `aws sso login --profile joebarbere-admin`; without a
   session the live leg reports SKIPPED, not OK.
5. Prove any new deny with a `--dry-run` call (`UnauthorizedOperation` = denied,
   `DryRunOperation` = allowed). A guard is tested through the path that runs it.

## Also live in the account (created by CLI 2026-09-26; import into Terraform later)

- Cost-allocation tags `Project` and `ManagedBy` — **active** (the split starts on that date).
- Budgets, alerts only, email at 50/80/100% actual and 100% forecast: `account-monthly-25`
  ($25/month, account-wide), `jansky-research-monthly` and `aws-ai-monthly` (per `Project`).
- Cost-anomaly monitor `services-anomaly-monitor`, daily email for anomalies ≥ $5.

No compute, storage or Terraform exists for this repo yet.
