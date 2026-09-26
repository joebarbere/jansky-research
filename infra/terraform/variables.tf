variable "region" {
  description = "Must be us-east-1: the aws-cloud seatbelt denies every other region."
  type        = string
  default     = "us-east-1"
}

variable "name_prefix" {
  description = "Prefix for resource names."
  type        = string
  default     = "jansky"
}

variable "gpu_instance_type" {
  description = "Default instance type in the GPU launch template. Must be on the aws-cloud seatbelt allowlist."
  type        = string
  default     = "g6.xlarge"
}

variable "gpu_ami_ssm_parameter" {
  description = "Public SSM parameter resolving to the NVIDIA-driver base AMI (driver + CUDA, no framework)."
  type        = string
  default     = "/aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-22.04/latest/ami-id"
}

variable "max_lifetime_minutes" {
  description = "Hard ceiling: user-data schedules a shutdown this long after boot, and shutdown means terminate. A forgotten g6.xlarge on-demand is ~$580/month; this caps it at a few dollars."
  type        = number
  default     = 120
}
