# One-shot compute: a launch template that jobs are started from with
# `aws ec2 run-instances --launch-template`. Nothing here bills until an instance runs.
#
# Access is SSM only (no SSH key, no inbound rules), the instance lives in the default VPC's
# public subnets (so no NAT gateway -- the seatbelt denies those anyway), and it terminates
# itself: shutdown behaviour is terminate, and user-data schedules a shutdown at boot.

data "aws_ssm_parameter" "gpu_ami" {
  name = var.gpu_ami_ssm_parameter
}

data "aws_vpc" "default" {
  default = true
}

data "aws_iam_policy_document" "ec2_trust" {
  statement {
    actions = ["sts:AssumeRole"]
    principals {
      type        = "Service"
      identifiers = ["ec2.amazonaws.com"]
    }
  }
}

resource "aws_iam_role" "compute" {
  name               = "${var.name_prefix}-compute"
  assume_role_policy = data.aws_iam_policy_document.ec2_trust.json
}

resource "aws_iam_role_policy_attachment" "ssm" {
  role       = aws_iam_role.compute.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_instance_profile" "compute" {
  name = "${var.name_prefix}-compute"
  role = aws_iam_role.compute.name
}

resource "aws_security_group" "compute" {
  name        = "${var.name_prefix}-compute"
  description = "No inbound. Outbound for package installs, the repo clone and SSM."
  vpc_id      = data.aws_vpc.default.id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_launch_template" "gpu" {
  name                                 = "${var.name_prefix}-gpu"
  image_id                             = data.aws_ssm_parameter.gpu_ami.value
  instance_type                        = var.gpu_instance_type
  instance_initiated_shutdown_behavior = "terminate"
  vpc_security_group_ids               = [aws_security_group.compute.id]
  update_default_version               = true

  iam_instance_profile {
    arn = aws_iam_instance_profile.compute.arn
  }

  metadata_options {
    http_tokens   = "required" # IMDSv2 only
    http_endpoint = "enabled"
  }

  block_device_mappings {
    device_name = "/dev/sda1"
    ebs {
      volume_size           = 75
      volume_type           = "gp3"
      delete_on_termination = true
      encrypted             = true
    }
  }

  user_data = base64encode(<<-EOT
    #!/bin/bash
    shutdown -h +${var.max_lifetime_minutes} "jansky max-lifetime reached"
  EOT
  )

  # default_tags do not reach instances launched from a template by the CLI; these do.
  dynamic "tag_specifications" {
    for_each = ["instance", "volume"]
    content {
      resource_type = tag_specifications.value
      tags = {
        Project   = "jansky-research"
        ManagedBy = "terraform"
        Purpose   = "research"
        Name      = "${var.name_prefix}-gpu-job"
      }
    }
  }
}
