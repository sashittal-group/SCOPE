import pandas as pd
import os

def main():
    df = pd.read_csv("data/williams/signatures_dataset/DLP/SNV/snv_counts.csv.gz", low_memory=False)
    df_cn_clones = pd.read_csv("data/williams/signatures_dataset/DLP/CNA/hscn_clones.csv.gz")

    patient_ids = df['patient'].unique().tolist()

    for patient_id in patient_ids:
        print(patient_id)
        os.makedirs(f"data/williams/scratch/{patient_id}", exist_ok=True)
        df_snv_counts = df[df['patient'] == patient_id]
        df_snv_counts.to_csv(f"data/williams/scratch/{patient_id}/snv_counts.csv", index=False)

        df_cn_clones_patient = df_cn_clones[df_cn_clones['patient'] == patient_id]
        df_cn_clones_patient.to_csv(f"data/williams/scratch/{patient_id}/hscn_clones.csv", index=False)
    
    with open("data/williams/scratch/separate.txt", 'w') as f:
        print("Done", file=f)


if __name__ == "__main__":
    main()