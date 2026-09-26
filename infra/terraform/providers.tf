# jansky-research's own AWS resources (plan 96). Account-wide guardrails -- the seatbelt,
# budgets, cost-allocation tags -- belong to the sibling aws-cloud repo, never here.
provider "aws" {
  region = var.region

  default_tags {
    tags = {
      Project   = "jansky-research"
      ManagedBy = "terraform"
      Purpose   = "research"
    }
  }
}
