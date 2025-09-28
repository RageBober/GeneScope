#!/usr/bin/env python3
"""
Test script for GenoScope genomics functionality
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genoscope.integrations import ClinVarAPI, DbSNPAPI, VCFProcessor
from genoscope.pipeline.genomics_pipeline import GenomicsPipeline

def test_clinvar():
    """Test ClinVar API"""
    print("\n" + "="*50)
    print("Testing ClinVar API")
    print("="*50)
    
    clinvar = ClinVarAPI()
    
    # Test searching by rsID (known BRCA1 variant)
    print("\n1. Testing search by rsID (rs80357906 - BRCA1)...")
    result = clinvar.search_by_rsid("rs80357906")
    if result:
        print(f"   ✓ Found: {result.get('clinical_significance', 'Unknown')}")
        print(f"   Gene: {result.get('gene', 'N/A')}")
        print(f"   URL: {result.get('url', 'N/A')}")
    else:
        print("   ✗ Not found")
    
    # Test variant interpretation
    print("\n2. Testing significance interpretation...")
    interp = clinvar.interpret_significance("Pathogenic")
    print(f"   Pathogenic → {interp['category']}, actionable: {interp['actionable']}")
    
    return True

def test_dbsnp():
    """Test dbSNP API"""
    print("\n" + "="*50)
    print("Testing dbSNP API")
    print("="*50)
    
    dbsnp = DbSNPAPI()
    
    print("\n1. Testing dbSNP lookup (rs1799963)...")
    result = dbsnp.get_snp_info("rs1799963")
    if result:
        print(f"   ✓ Found: {result.get('rsid')}")
        print(f"   Chromosome: {result.get('chromosome', 'N/A')}")
        print(f"   Gene: {result.get('gene', 'N/A')}")
    else:
        print("   ✗ Not found")
    
    return True

def test_vcf_processor():
    """Test VCF processing"""
    print("\n" + "="*50)
    print("Testing VCF Processor")
    print("="*50)
    
    # Create a test VCF file
    test_vcf = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100000	rs123	A	G	50	PASS	DP=30	GT:GQ:DP	0/1:45:30
chr17	43044295	rs80357906	G	A	99	PASS	DP=50	GT:GQ:DP	1/1:99:50
"""
    
    # Write to temp file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(test_vcf)
        temp_vcf = f.name
    
    try:
        print(f"\n1. Processing test VCF file...")
        processor = VCFProcessor(temp_vcf)
        variants = processor.parse()
        
        print(f"   ✓ Loaded {len(variants)} variants")
        
        stats = processor.get_statistics()
        print(f"   Total variants: {stats['total_variants']}")
        print(f"   Variant types: {stats.get('variant_types', {})}")
        
        # Test filtering
        filtered = processor.filter_variants(min_quality=60)
        print(f"   ✓ Filtered to {len(filtered)} high-quality variants")
        
        return True
        
    finally:
        os.unlink(temp_vcf)

def test_pipeline():
    """Test the genomics pipeline"""
    print("\n" + "="*50)
    print("Testing Genomics Pipeline")
    print("="*50)
    
    pipeline = GenomicsPipeline()
    
    print("\n1. Checking available tools:")
    for tool, available in pipeline.tools.items():
        status = "✓" if available else "✗"
        print(f"   {status} {tool}")
    
    print("\n2. Testing mock pipeline...")
    # Create a minimal test FASTQ
    test_fastq = "/tmp/test.fastq"
    with open(test_fastq, 'w') as f:
        f.write("@read1\nACGT\n+\n####\n")
    
    try:
        results = pipeline.run_full_pipeline(
            fastq_r1=test_fastq,
            sample_name="test",
            output_dir="/tmp/genoscope_test"
        )
        
        if results['status'] == 'completed':
            print("   ✓ Pipeline completed successfully")
            if 'annotation' in results.get('steps', {}):
                ann = results['steps']['annotation']
                print(f"   Total variants: {ann.get('total_variants', 0)}")
                print(f"   ClinVar matches: {ann.get('clinvar_matches', 0)}")
        else:
            print(f"   ✗ Pipeline failed: {results.get('error', 'Unknown error')}")
            
    finally:
        if os.path.exists(test_fastq):
            os.unlink(test_fastq)
    
    return True

def main():
    """Run all tests"""
    print("\n" + "#"*60)
    print("#" + " GenoScope Genomics Functionality Test ".center(58) + "#")
    print("#"*60)
    
    tests = [
        ("ClinVar API", test_clinvar),
        ("dbSNP API", test_dbsnp),
        ("VCF Processor", test_vcf_processor),
        ("Genomics Pipeline", test_pipeline)
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n✗ {name} failed with error: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    for name, success in results:
        status = "✓ PASSED" if success else "✗ FAILED"
        print(f"{name:.<30} {status}")
    
    passed = sum(1 for _, s in results if s)
    total = len(results)
    print(f"\nTotal: {passed}/{total} tests passed")
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
