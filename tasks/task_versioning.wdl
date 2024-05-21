version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
  }
  meta {
    volatile: true
  }
  command {
    WF_Version="v1.0.3 (PHB v1.2.1)"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$WF_Version" > WF_VERSION
  }
  output {
    String date = read_string("TODAY")
    String wf_version = read_string("WF_VERSION")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: docker
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}

