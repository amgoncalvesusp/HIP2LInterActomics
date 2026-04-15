import os

def processar_docking(pasta_origem, last_pa):
    pasta_prot = 'proteinas_pdb'
    pasta_lig = 'ligantes_mol2'
    
    for p in [pasta_prot, pasta_lig]:
        if not os.path.exists(p): os.makedirs(p)

    arquivos = [f for f in os.listdir(pasta_origem) if f.endswith('.mol2')]

    for nome_arq in arquivos:
        caminho = os.path.join(pasta_origem, nome_arq)
        nome_base = os.path.splitext(nome_arq)[0]
        
        with open(caminho, 'r') as f:
            linhas = f.readlines()

        proteina_atoms = []
        
        # Usaremos listas para guardar os DADOS brutos antes de formatar
        ligante_atoms_data = [] 
        ligante_bonds_data = []
        
        secao_atual = None
        
        # Se last_pa é o último da proteína (ex: 4068), o ligante começa em last_pa + 1
        primeiro_id_ligante = last_pa + 1  

        for linha in linhas:
            if "@<TRIPOS>ATOM" in linha:
                secao_atual = "ATOM"
                continue
            elif "@<TRIPOS>BOND" in linha:
                secao_atual = "BOND"
                continue
            elif "@<TRIPOS>SUBSTRUCTURE" in linha:
                secao_atual = "SUB"
                continue

            if secao_atual == "ATOM":
                partes = linha.split()
                if len(partes) < 2: continue
                idx = int(partes[0])
                
                if idx <= last_pa:
                    # Formata para PDB
                    try:
                        pdb_line = f"ATOM  {idx:>5} {partes[1]:^4} {partes[7][:3]:>3} A{partes[7][3:]:>4}    {float(partes[2]):8.3f}{float(partes[3]):8.3f}{float(partes[4]):8.3f}  1.00  0.00\n"
                        proteina_atoms.append(pdb_line)
                    except: pass
                else:
                    # É ligante - Calculamos o novo ID
                    novo_idx = idx - last_pa
                    
                    # Salvamos o "resto" da linha intacto para não perder o alinhamento original das coordenadas
                    idx_str = partes[0]
                    idx_end_pos = linha.find(idx_str) + len(idx_str)
                    resto_da_linha = linha[idx_end_pos:] 
                    
                    ligante_atoms_data.append((novo_idx, resto_da_linha))

            elif secao_atual == "BOND":
                partes = linha.split()
                if len(partes) < 3: continue
                at1, at2 = int(partes[1]), int(partes[2])
                
                # Só pega a ligação se ambos os átomos forem do ligante
                if at1 >= primeiro_id_ligante and at2 >= primeiro_id_ligante:
                    novo_at1 = at1 - last_pa
                    novo_at2 = at2 - last_pa
                    b_idx = len(ligante_bonds_data) + 1
                    
                    ligante_bonds_data.append((b_idx, novo_at1, novo_at2, partes[3]))

        # --- ETAPA DE FORMATAÇÃO DO LIGANTE ---
        num_atoms = len(ligante_atoms_data)
        num_bonds = len(ligante_bonds_data)
        
        # Calcula os espaços ocupados (EOIA e EOLA)
        EOIA = len(str(num_atoms))
        EOLA = len(str(num_bonds))

        ligante_atoms_formatados = []
        for novo_idx, resto_da_linha in ligante_atoms_data:
            # 2 espaços vazios + índice alinhado à direita pelo tamanho do EOIA + resto da linha
            linha_formatada = f"  {str(novo_idx).rjust(EOIA)}{resto_da_linha}"
            ligante_atoms_formatados.append(linha_formatada)

        ligante_bonds_formatados = []
        for b_idx, at1, at2, tipo in ligante_bonds_data:
            # 2 espaços vazios + index_bond(EOLA) + espaço + at1(EOIA) + espaço + at2(EOIA) + tipo
            linha_formatada = f"  {str(b_idx).rjust(EOLA)} {str(at1).rjust(EOIA)} {str(at2).rjust(EOIA)} {tipo:>4}\n"
            ligante_bonds_formatados.append(linha_formatada)

        # --- SALVAR PROTEÍNA PDB ---
        if proteina_atoms:
            with open(os.path.join(pasta_prot, f"{nome_base}.pdb"), 'w') as f_p:
                f_p.write("REMARK   Separated Protein\n")
                f_p.writelines(proteina_atoms)
                f_p.write("END\n")

        # --- SALVAR LIGANTE MOL2 ---
        if ligante_atoms_formatados:
            with open(os.path.join(pasta_lig, f"{nome_base}_ligand.mol2"), 'w') as f_l:
                f_l.write("@<TRIPOS>MOLECULE\n")
                f_l.write(f"{nome_base}_LIG\n")
                f_l.write(f"{num_atoms} {num_bonds} 1 0 0\n")
                f_l.write("SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
                f_l.writelines(ligante_atoms_formatados)
                f_l.write("@<TRIPOS>BOND\n")
                f_l.writelines(ligante_bonds_formatados)

        print(f"Sucesso: {nome_arq} | Átomos: {num_atoms} | Bonds: {num_bonds}")

last_pa = int(input('Escreva o número do último átomo da proteína: '))
processar_docking('mol2', last_pa)