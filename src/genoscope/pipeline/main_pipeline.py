    # Find input files
    input_dir = Path(config.input_dir)
    fastq_files = list(input_dir.glob("*.fastq*"))
    
    if len(fastq_files) == 0:
        raise ValueError(f"No FASTQ files found in {input_dir}")
    
    # Identify R1 and R2
    r1_files = [f for f in fastq_files if "_R1" in f.name or "_1." in f.name]
    r2_files = [f for f in fastq_files if "_R2" in f.name or "_2." in f.name]
    
    if r1_files:
        fastq_r1 = str(r1_files[0])
        fastq_r2 = str(r2_files[0]) if r2_files else None
    else:
        fastq_r1 = str(fastq_files[0])
        fastq_r2 = str(fastq_files[1]) if len(fastq_files) > 1 else None
    
    # Run pipeline
    return pipeline.run(fastq_r1, fastq_r2)
