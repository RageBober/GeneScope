"""
🔬 BioForge Workflow Examples

This module demonstrates how to use bioinformatics tools in integrated workflows.

Each example shows a real-world use case combining multiple tools.
"""

from pathlib import Path
from bioinformatics_tools import (
    MetagenomicsPipeline,
    ProteinAnalysisPipeline,
    MutationEffectPipeline,
    MetaGraphDataSource,
)


# ═══════════════════════════════════════════════════════════════
#  Example 1: Complete Metagenomics Workflow
# ═══════════════════════════════════════════════════════════════

def example_metagenomics_analysis():
    """
    Complete metagenomics pipeline:
    1. Classify organisms with Kraken2
    2. Assemble genomes with MEGAHIT
    3. Normalize taxonomy with GTDB-Tk

    Use Case: Microbiome study from human gut samples
    """
    print("=" * 60)
    print("Example 1: Metagenomics Analysis")
    print("=" * 60)

    # Initialize pipeline
    pipeline = MetagenomicsPipeline(
        kraken2_db="/data/kraken2_standard",
        gtdbtk_db="/data/gtdbtk_r214",
        assembler="megahit",  # Fast and memory-efficient
        use_docker=True
    )

    # Run complete workflow
    result = pipeline.run(
        fastq_r1="data/gut_sample_R1.fastq.gz",
        fastq_r2="data/gut_sample_R2.fastq.gz",
        output_dir="results/metagenomics"
    )

    # Print results
    print(f"\nPipeline Status: {result.status}")
    print(f"Steps Completed: {', '.join(result.steps_completed)}")

    if "kraken2" in result.results:
        print(f"\n📊 Taxonomic Classification:")
        print(f"  Classified: {result.results['kraken2']['classified_percentage']:.1f}%")
        print(f"  Total reads: {result.results['kraken2']['total_reads']:,}")

        print(f"\n🦠 Top 5 Species:")
        for i, taxon in enumerate(result.results.get('kraken2_top_species', [])[:5], 1):
            print(f"  {i}. {taxon['name']}: {taxon['percentage']:.2f}%")

    if "assembly" in result.results:
        print(f"\n🧬 Genome Assembly:")
        print(f"  Contigs: {result.results['assembly']['num_contigs']:,}")
        print(f"  Total length: {result.results['assembly']['total_length']:,} bp")
        print(f"  N50: {result.results['assembly']['n50']:,} bp")
        print(f"  Longest contig: {result.results['assembly']['longest_contig']:,} bp")

    # Save results
    output_path = Path("results/metagenomics/pipeline_result.json")
    result.to_json(output_path)
    print(f"\n✅ Results saved to: {output_path}")

    return result


# ═══════════════════════════════════════════════════════════════
#  Example 2: Protein Function Prediction
# ═══════════════════════════════════════════════════════════════

def example_protein_analysis():
    """
    Protein property prediction workflow:
    1. Analyze amino acid sequences with XtriMoPGLM
    2. Predict function, structure, and properties
    3. Prepare data for evolutionary simulations

    Use Case: Characterize novel proteins from metagenomic assembly
    """
    print("\n" + "=" * 60)
    print("Example 2: Protein Analysis")
    print("=" * 60)

    # Initialize pipeline
    pipeline = ProteinAnalysisPipeline(
        device="cuda",  # Use GPU for faster inference
        batch_size=64,
        use_docker=False  # Run locally for GPU access
    )

    # Run analysis
    result = pipeline.run(
        protein_fasta="data/novel_proteins.fasta",
        output_dir="results/protein_analysis"
    )

    print(f"\nPipeline Status: {result.status}")

    if "predictions" in result.results:
        predictions = result.results['predictions']
        print(f"\n🧪 Protein Predictions:")
        print(f"  Total proteins analyzed: {len(predictions)}")
        print(f"  Predictions include: function, structure, properties")

    return result


# ═══════════════════════════════════════════════════════════════
#  Example 3: Mutation Effect Prediction for EvoSim
# ═══════════════════════════════════════════════════════════════

def example_mutation_effects():
    """
    Mutation effect prediction workflow:
    1. Predict phenotypic effects of DNA mutations with Enformer
    2. Forecast gene expression changes
    3. Use predictions in EvoSim for evolutionary modeling

    Use Case: Study adaptive mutations in bacterial evolution
    """
    print("\n" + "=" * 60)
    print("Example 3: Mutation Effect Prediction")
    print("=" * 60)

    # Initialize pipeline
    pipeline = MutationEffectPipeline(
        model_name="enformer",  # or "deepsea"
        device="cuda",
        use_docker=False
    )

    # Run prediction
    result = pipeline.run(
        sequence_fasta="data/candidate_mutations.fasta",
        output_dir="results/mutation_effects"
    )

    print(f"\nPipeline Status: {result.status}")

    if "predictions" in result.results:
        predictions = result.results['predictions']
        print(f"\n🧬 Mutation Effects:")
        print(f"  Mutations analyzed: {len(predictions)}")
        print(f"  Effect predictions available for EvoSim integration")

    return result


# ═══════════════════════════════════════════════════════════════
#  Example 4: MetaGraph Global Sequence Search
# ═══════════════════════════════════════════════════════════════

def example_metagraph_search():
    """
    MetaGraph sequence discovery workflow:
    1. Search for sequences across global databases
    2. Find homologous genes in public datasets
    3. Integrate found sequences into GeneScope analysis

    Use Case: Discover related genes in SRA (Sequence Read Archive)
    """
    print("\n" + "=" * 60)
    print("Example 4: MetaGraph Sequence Search")
    print("=" * 60)

    # Initialize data source
    data_source = MetaGraphDataSource(
        index_dir="/data/metagraph_sra_index",
        use_docker=True
    )

    # Search for sequences
    result = data_source.search_sequences(
        query_fasta="data/query_genes.fasta",
        output_dir="results/metagraph_search",
        discovery_fraction=0.8
    )

    print(f"\nSearch Status: {result.status}")

    if "num_matches" in result.results:
        print(f"\n🔍 Search Results:")
        print(f"  Matches found: {result.results['num_matches']}")
        print(f"  Matched datasets available for download")

    return result


# ═══════════════════════════════════════════════════════════════
#  Example 5: Combined Workflow - Metagenomics + Protein Analysis
# ═══════════════════════════════════════════════════════════════

def example_combined_workflow():
    """
    Multi-stage workflow combining several pipelines:
    1. Run metagenomics analysis to get assembled genomes
    2. Extract protein sequences from assemblies
    3. Predict protein functions with XtriMoPGLM
    4. Analyze mutations with Enformer

    Use Case: Complete characterization of novel microbial community
    """
    print("\n" + "=" * 60)
    print("Example 5: Combined Metagenomics + Protein Analysis")
    print("=" * 60)

    # Stage 1: Metagenomics
    print("\n📍 Stage 1/3: Metagenomics Analysis")
    metagenomics = MetagenomicsPipeline(
        kraken2_db="/data/kraken2_standard",
        assembler="megahit"
    )

    mg_result = metagenomics.run(
        fastq_r1="data/sample_R1.fastq.gz",
        fastq_r2="data/sample_R2.fastq.gz",
        output_dir="results/combined_workflow/metagenomics"
    )

    print(f"✅ Metagenomics: {mg_result.status}")

    # Stage 2: Protein Prediction
    # (In real workflow, extract proteins from contigs first using tools like Prodigal)
    print("\n📍 Stage 2/3: Protein Analysis")
    protein_pipeline = ProteinAnalysisPipeline(device="cuda")

    protein_result = protein_pipeline.run(
        protein_fasta="results/combined_workflow/metagenomics/predicted_proteins.fasta",
        output_dir="results/combined_workflow/proteins"
    )

    print(f"✅ Protein Analysis: {protein_result.status}")

    # Stage 3: Mutation Effects
    print("\n📍 Stage 3/3: Mutation Effect Prediction")
    mutation_pipeline = MutationEffectPipeline(model_name="enformer")

    mutation_result = mutation_pipeline.run(
        sequence_fasta="results/combined_workflow/metagenomics/variant_regions.fasta",
        output_dir="results/combined_workflow/mutations"
    )

    print(f"✅ Mutation Effects: {mutation_result.status}")

    # Summary
    print("\n" + "=" * 60)
    print("📊 Combined Workflow Summary")
    print("=" * 60)

    if "assembly" in mg_result.results:
        print(f"\n🧬 Assembly: {mg_result.results['assembly']['num_contigs']} contigs")

    if "kraken2_top_species" in mg_result.results:
        print(f"\n🦠 Top organism: {mg_result.results['kraken2_top_species'][0]['name']}")

    print(f"\n✅ All stages completed successfully!")

    return {
        "metagenomics": mg_result,
        "proteins": protein_result,
        "mutations": mutation_result,
    }


# ═══════════════════════════════════════════════════════════════
#  Main Runner
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    """
    Run all examples (requires databases and sample data).

    To run specific examples:
        python bioforge_workflows.py
    """
    import sys

    print("""
╔══════════════════════════════════════════════════════════════╗
║            BioForge Workflow Examples                        ║
║                                                              ║
║  Demonstrating integrated bioinformatics pipelines          ║
╚══════════════════════════════════════════════════════════════╝
    """)

    # Check if databases are available
    print("\n⚠️  Note: These examples require:")
    print("  - Kraken2 database (~50GB)")
    print("  - GTDB-Tk database (~80GB)")
    print("  - MetaGraph index (varies)")
    print("  - Sample FASTQ files")
    print("  - GPU for ML models (optional)")

    choice = input("\nRun examples in DRY RUN mode (no actual execution)? [Y/n]: ")

    if choice.lower() in ['n', 'no']:
        print("\n⚠️  Running with actual tool execution...")
        print("Make sure all databases and sample data are available!\n")
    else:
        print("\n✅ Running in DRY RUN mode (showing workflow structure only)\n")
        # Set dry_run mode globally here if needed

    try:
        # Run examples
        print("\n" + "🚀 Starting workflow examples..." + "\n")

        # Example 1: Metagenomics
        # example_metagenomics_analysis()

        # Example 2: Protein Analysis
        # example_protein_analysis()

        # Example 3: Mutation Effects
        # example_mutation_effects()

        # Example 4: MetaGraph Search
        # example_metagraph_search()

        # Example 5: Combined Workflow
        # example_combined_workflow()

        print("\n" + "=" * 60)
        print("📚 Example workflows demonstrated successfully!")
        print("=" * 60)
        print("\nTo run actual workflows:")
        print("  1. Download required databases")
        print("  2. Prepare sample data")
        print("  3. Uncomment example calls in main()")
        print("  4. python bioforge_workflows.py")

    except Exception as exc:
        print(f"\n❌ Error: {exc}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
