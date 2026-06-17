import hipplinteractomics_terminal
import json
import os
import sys
import subprocess
from pathlib import Path

FOLDERS = [
"a3b4_anh_asp","a3b4_anh_chemplp","a3b4_anh_chemscore","a3b4_anh_goldscore"
,"a3b4_hyd_asp","a3b4_hyd_chemplp","a3b4_hyd_chemscore","a3b4_hyd_goldscore"
,"a4b2_anh_asp","a4b2_anh_chemplp","a4b2_anh_chemscore","a4b2_anh_goldscore"
,"a4b2_hyd_asp","a4b2_hyd_chemplp","a4b2_hyd_chemscore","a4b2_hyd_goldscore"
,"a7o_anh_asp","a7o_anh_chemplp","a7o_anh_chemscore","a7o_anh_goldscore"
,"a7o_hyd_asp","a7o_hyd_chemplp","a7o_hyd_chemscore","a7o_hyd_goldscore"
,"a7t_anh_asp","a7t_anh_chemplp","a7t_anh_chemscore","a7t_anh_goldscore"
]

if  __name__ == "__main__":
    ### setting params

    fp_length = [4096,2048,1024] # test :[1024]#
    fp_rd_lv =  [(2,10),(3,5),(6,2)] # somando 10 A esfera maxima  # test :[(3,5)]#
    fp_format =  ['bin','cnt'] # test: ['bin']#
    app_path= "/home/laqmedsom/softwares_laqmedsomm/HIP2LInterActomics_GUI/HIP2LInterActomics_GUI-terminal-fp-trajectory-updates-20260521"
    work_path = "/home/laqmedsom/Documentos/Daniel/PhDProject_correction/6-Calc_IFPS_correct/inputs"
    path_ =  [f'{work_path}/{i}' for i in FOLDERS] # test: [f'{work_path}/a4b2_hyd_chemplp'] #
    
    activity_rep = {"a3b4":"pKi_value","a4b2":"pKi_value","a7o":"pKi_value","a7t":"pEC50_value"}
    ref_activity_files = {"a3b4":"2-a3b4_LigandsToGeometryOptimization.tsv","a4b2":"2-a4b2_SH-EP1_LigandsToGeometryOptimization.tsv",
                          "a7o":"2-a7o_LigandsToGeometryOptimization.tsv","a7t":"2-a7t_ec50_LigandsToGeometryOptimization.tsv"}
    
    ## program ########
    os.system(f"python {app_path}/hipplinteractomics_terminal.py --write-template {work_path}/hipplinteractomics_terminal_conf_example.txt")
    
    nproc_per_job = 6
    max_parallel = 2
    
    in_dict={}
    commands = []
    with open(f'{work_path}/hipplinteractomics_terminal_conf_example.txt', 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line :
                continue
            if line.count(':') >= 1: 
                key,value = line.split(':')
            else: 
                continue    
            
            in_dict[key.strip().strip('"')] = value.strip().strip('"').strip('",')
            
    running = []        
    
    for path in path_: 
        for length in fp_length:
            for rd_l in fp_rd_lv:
                for form in fp_format:
                    folder_name = path.split('/')[-1]
                    site= folder_name.split("_")[0]
                    
                    bcode = f"B{length:04d}"
                    lrcode = f"R{rd_l[0]}_L{rd_l[1]}"
                    work_dir = f"{path}/{folder_name}_{bcode}_{form}_{lrcode}"  
                    
                    os.makedirs(work_dir, exist_ok=True)
                                 
                    in_dict["protein_file"]        = f"{path}/proteinas_pdb"
                    in_dict["ligand_file"]         = f"{path}/ligantes_mol2"
                    in_dict["workdir"]             = work_dir
                    in_dict["ifp_levels"]          = rd_l[1]
                    in_dict["ifp_radius"]          = rd_l[0]
                    in_dict["ifp_length"]          = length
                    in_dict["ifp_bit"]             = True if form == 'bin' else False
                    in_dict["nproc"]               = nproc_per_job
                    in_dict["fp_labels_csv"]       = f'{work_path}/{ref_activity_files[site]}'
                    in_dict["fp_labels_id_column"] = 'molecule_chembl_id'
                    in_dict["fp_labels_column"]    = activity_rep[site]

                    print(in_dict)
                    conf_name = f"hipplinteractomics_terminal_conf_{folder_name}_{bcode}_{form}_{lrcode}.json"
                    conf_path = f"{work_dir}/{conf_name}"
                    
                    with open(conf_path, "w", encoding = "utf-8") as file:
                        json.dump(in_dict, file)
                        
                    log_path = f"{work_dir}/run_{folder_name}_{bcode}_{form}_{lrcode}.log"    
                    
                    cmd = [
                        sys.executable
                        ,f"{app_path}/hipplinteractomics_terminal.py"
                        ,conf_path
                    ]
                    
                    log_file = open(log_path, "w", encoding="utf-8")

                    print("Starting:", folder_name, bcode, form, lrcode)

                    p = subprocess.Popen(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT
                    )

                    running.append((p, log_file, log_path))

                    if len(running) >= max_parallel:
                        p_old, log_old, log_old_path = running.pop(0)
                        return_code = p_old.wait()
                        log_old.close()

                        print("Finished:", log_old_path, "return code:", return_code)

    for p, log_file, log_path in running:
        return_code = p.wait()
        log_file.close()
        
        print("Finished:", log_path, "return code:", return_code)
          
        print("All jobs finished.")