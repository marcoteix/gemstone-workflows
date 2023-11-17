version 1.0

task taxonomy_qc_check {
    input {
        String samplename
        String gambit_taxonomy
        Float checkm2_contamination
        String? lab_determined_genus
        Float contamination_threshold = 0.01
        Int disk_size = 20
    }
    command <<<
        python3 <<CODE
        import os 

        lab_determined_genus = "~{lab_determined_genus}"
        checkm2_contamination = ~{checkm2_contamination}
        gambit_taxonomy = "~{gambit_taxonomy}"
        contamination_threshold = ~{contamination_threshold}
        samplename = "~{samplename}"

        # Take the genus from the taxonomy predictions
        valid_lab_genus = True
        if type(lab_determined_genus) == str:
            lab_determined_genus = lab_determined_genus.replace('_', ' ').split(' ')[0].capitalize()
            # Sometimes the lab determined genus has 2 genera (usually Escherichia/Shigella)
            # separated by a /. In these cases, include both genera in the comparison
            lab_determined_genus = [x.capitalize() for x in lab_determined_genus.split('/')]
            # Do not fail QC if the lab determined genus is unknown. Check if lab predicted genus matched the GAMBIT prediction
            if lab_determined_genus is None: valid_lab_genus = False
            if any([x in ['', 'Unknown'] for x in lab_determined_genus]): valid_lab_genus = False
        else: valid_lab_genus = False

        gambit_taxonomy = gambit_taxonomy.replace('_', ' ').split(' ')[0].capitalize()

        lab_genus_qc_check = 'QC_ALERT'
        contamination_qc_check = 'QC_ALERT'


        if valid_lab_genus:
            if gambit_taxonomy in lab_determined_genus:
                print(f'DEBUG: Found matching genera between the lab and GAMBIT predictions ({" or ".join(lab_determined_genus)} and \
        {gambit_taxonomy}, respectively). Setting QC flag to QC_PASS.')
                lab_genus_qc_check = 'QC_PASS'
            else:
                print(f'DEBUG: Found no matching genera between the lab and GAMBIT predictions ({" or ".join(lab_determined_genus)} and \
        {gambit_taxonomy}, respectively). Setting QC flag to QC_ALERT.')
                lab_genus_qc_check = 'QC_ALERT'
        else:
            print(f'DEBUG: Got no lab predicted genera (or prediction is NaN or unknown). Setting the QC flag to QC_PASS.')
            lab_genus_qc_check = 'QC_PASS'

        # Check if the checkM2 contamination ratio is less than the threshold
        if checkm2_contamination <= contamination_threshold:
            print(f'DEBUG: The checkM2 contamination value ({checkm2_contamination}) is less than or equal to the set contamination limit \
        ({contamination_threshold}). Setting QC flag to QC_PASS.')
            contamination_qc_check = 'QC_PASS'
        else:
            print(f'DEBUG: The checkM2 contamination value ({checkm2_contamination}) is greater than the set contamination limit \
        ({contamination_threshold}). Setting QC flag to QC_ALERT.')
            contamination_qc_check = 'QC_ALERT'

        if lab_genus_qc_check == 'QC_ALERT' or contamination_qc_check == 'QC_ALERT': qc_flag = 'QC_ALERT'
        else: qc_flag = 'QC_PASS'

        print(f'DEBUG: Global QC: {qc_flag} | Lab genus QC: {lab_genus_qc_check} | Contamination QC: {contamination_qc_check}')

        # Make dir if it does not exist
        if not os.path.exists(samplename): os.mkdir(samplename)

        with open(f"{samplename}/QC_FLAG.txt", 'w') as f:
            f.write(qc_flag)
        with open(f"{samplename}/QC_report.txt", 'w') as f:
            f.write(f'Global QC: {qc_flag}\nLab genus QC: {lab_genus_qc_check}\nContamination QC: {contamination_qc_check}')

        CODE
    >>>
    output {
        String qc_check = read_string("~{samplename}/QC_FLAG.txt")
        File qc_report = "~{samplename}/QC_report.txt"
    }
    runtime {
        docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
        preemptible: 0
    }
}