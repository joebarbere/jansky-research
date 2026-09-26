output "gpu_launch_template" {
  description = "Launch a one-shot GPU job with: aws ec2 run-instances --launch-template LaunchTemplateName=<this>"
  value       = aws_launch_template.gpu.name
}

output "gpu_ami" {
  description = "The AMI the template resolved to at apply time."
  value       = nonsensitive(data.aws_ssm_parameter.gpu_ami.value)
}
