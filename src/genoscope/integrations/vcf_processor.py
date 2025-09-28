"""
VCF File Processing Module
Handles parsing and processing of VCF (Variant Call Format) files
"""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Union
from pathlib import Path
import gzip
import re


@dataclass
class Variant:
    """Represents a single genetic variant from VCF file"""
    chromosome: str
    position: int
    rsid: Optional[str]
    ref: str
    alt: str
    quality: Optional[float]
    filter_status: str
    info: Dict[str, Any]
    format_fields: Optional[str] = None
    genotype: Optional[str] = None
    
    def get_variant_type(self) -> str:
        """Determine the type of variant (SNP, insertion, deletion, etc.)"""
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNP"
        elif len(self.ref) > len(self.alt):
            return "deletion" 
        elif len(self.ref) < len(self.alt):
            return "insertion"
        else:
            return "complex"
    
    def is_transition(self) -> bool:
        """Check if SNP is a transition (A<->G or C<->T)"""
        if self.get_variant_type() != "SNP":
            return False
        transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        return (self.ref, self.alt) in transitions


class VCFProcessor:
    """Process and analyze VCF files"""
    
    def __init__(self, vcf_path: Union[str, Path]):
        self.vcf_path = Path(vcf_path)
        self.variants: List[Variant] = []
        self.header_lines: List[str] = []
        self.sample_names: List[str] = []
        
        if not self.vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    def _open_file(self):
        """Open VCF file, handling both plain and gzipped files"""
        if self.vcf_path.suffix == '.gz':
            return gzip.open(self.vcf_path, 'rt')
        else:
            return open(self.vcf_path, 'r')
    
    def parse(self, limit: Optional[int] = None) -> List[Variant]:
        """Parse VCF file and return list of variants"""
        variants = []
        
        with self._open_file() as f:
            # Parse header
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    self.header_lines.append(line)
                elif line.startswith('#CHROM'):
                    # Column header line
                    columns = line.split('\t')
                    if len(columns) > 9:
                        self.sample_names = columns[9:]  # Sample names start from column 9
                    break
            
            # Parse variant lines
            for i, line in enumerate(f):
                if limit and i >= limit:
                    break
                    
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                variant = self._parse_variant_line(line)
                if variant:
                    variants.append(variant)
        
        self.variants = variants
        return variants
    
    def _parse_variant_line(self, line: str) -> Optional[Variant]:
        """Parse a single VCF variant line"""
        try:
            fields = line.split('\t')
            if len(fields) < 8:
                return None
            
            # Parse basic fields
            chrom = fields[0]
            pos = int(fields[1])
            rsid = fields[2] if fields[2] != '.' else None
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5]) if fields[5] != '.' else None
            filter_status = fields[6]
            info_str = fields[7]
            
            # Parse INFO field
            info = self._parse_info_field(info_str)
            
            # Parse genotype if available
            format_field = fields[8] if len(fields) > 8 else None
            genotype = fields[9] if len(fields) > 9 else None
            
            return Variant(
                chromosome=chrom,
                position=pos,
                rsid=rsid,
                ref=ref,
                alt=alt,
                quality=qual,
                filter_status=filter_status,
                info=info,
                format_fields=format_field,
                genotype=genotype
            )
            
        except (ValueError, IndexError) as e:
            # Skip malformed lines
            return None
    
    def _parse_info_field(self, info_str: str) -> Dict[str, Any]:
        """Parse the INFO field into a dictionary"""
        info = {}
        if info_str == '.':
            return info
        
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                # Try to convert to number if possible
                try:
                    if '.' in value:
                        info[key] = float(value)
                    else:
                        info[key] = int(value)
                except ValueError:
                    info[key] = value
            else:
                # Flag without value
                info[item] = True
        
        return info
    
    def filter_variants(self, min_quality: Optional[float] = None, 
                       pass_only: bool = False,
                       chromosomes: Optional[List[str]] = None) -> List[Variant]:
        """Filter variants based on criteria"""
        filtered = self.variants.copy()
        
        if min_quality is not None:
            filtered = [v for v in filtered if v.quality and v.quality >= min_quality]
        
        if pass_only:
            filtered = [v for v in filtered if v.filter_status == 'PASS']
        
        if chromosomes:
            filtered = [v for v in filtered if v.chromosome in chromosomes]
        
        return filtered
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get statistics about the variants"""
        if not self.variants:
            return {"total_variants": 0}
        
        variant_types = {}
        filter_status = {}
        quality_scores = []
        chromosomes = set()
        
        for variant in self.variants:
            # Count variant types
            vtype = variant.get_variant_type()
            variant_types[vtype] = variant_types.get(vtype, 0) + 1
            
            # Count filter status
            filter_status[variant.filter_status] = filter_status.get(variant.filter_status, 0) + 1
            
            # Collect quality scores
            if variant.quality:
                quality_scores.append(variant.quality)
            
            # Track chromosomes
            chromosomes.add(variant.chromosome)
        
        stats = {
            "total_variants": len(self.variants),
            "variant_types": variant_types,
            "filter_status": filter_status,
            "chromosomes": sorted(list(chromosomes)),
            "sample_count": len(self.sample_names)
        }
        
        if quality_scores:
            stats["quality_stats"] = {
                "min": min(quality_scores),
                "max": max(quality_scores),
                "mean": sum(quality_scores) / len(quality_scores)
            }
        
        return stats
    
    def export_to_dict(self) -> List[Dict[str, Any]]:
        """Export variants to list of dictionaries"""
        return [
            {
                "chromosome": v.chromosome,
                "position": v.position,
                "rsid": v.rsid,
                "ref": v.ref,
                "alt": v.alt,
                "quality": v.quality,
                "filter": v.filter_status,
                "variant_type": v.get_variant_type(),
                "info": v.info
            }
            for v in self.variants
        ]


# Example usage and testing
def main():
    """Test the VCF processor with a sample file"""
    
    # Create a test VCF file
    test_vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100000	rs123	A	G	50.0	PASS	DP=30	GT:GQ:DP	0/1:45:30
chr17	43044295	rs80357906	G	A	99.0	PASS	DP=50	GT:GQ:DP	1/1:99:50
chr2	200000	.	AT	A	25.5	LowQual	DP=15	GT:GQ:DP	0/1:20:15
"""
    
    # Write to temp file
    import tempfile
    import os
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(test_vcf_content)
        temp_vcf_path = f.name
    
    try:
        print("Testing VCF Processor")
        print("=" * 50)
        
        # Process VCF
        processor = VCFProcessor(temp_vcf_path)
        variants = processor.parse()
        
        print(f"Loaded {len(variants)} variants")
        
        # Get statistics
        stats = processor.get_statistics()
        print(f"Statistics:")
        print(f"  Total variants: {stats['total_variants']}")
        print(f"  Variant types: {stats['variant_types']}")
        print(f"  Filter status: {stats['filter_status']}")
        print()
        
        # Filter variants
        filtered = processor.filter_variants(min_quality=30, pass_only=True)
        print(f"Filtered to {len(filtered)} high-quality variants")
        print()
        
        # Show first variant details
        if variants:
            var = variants[0]
            print(f"First variant details:")
            print(f"  Position: {var.chromosome}:{var.position}")
            print(f"  Change: {var.ref} → {var.alt}")
            print(f"  Type: {var.get_variant_type()}")
            print(f"  Quality: {var.quality}")
            print(f"  Genotype: {var.genotype}")
        
    finally:
        # Clean up temp file
        os.unlink(temp_vcf_path)


if __name__ == "__main__":
    main()