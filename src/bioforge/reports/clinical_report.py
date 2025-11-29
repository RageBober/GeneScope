"""
Clinical Report Generator for Genomics Analysis
Generates professional PDF/HTML reports for clinical and research use
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import json
from jinja2 import Environment, FileSystemLoader
import logging

logger = logging.getLogger(__name__)

class ClinicalReportGenerator:
    """
    Generates clinical-grade genomics reports
    """
    
    def __init__(self, template_dir: Optional[str] = None):
        """
        Initialize report generator
        
        Args:
            template_dir: Directory containing report templates
        """
        if template_dir:
            self.template_dir = Path(template_dir)
        else:
            self.template_dir = Path(__file__).parent / "templates"
        
        # Create template directories if they don't exist
        self.template_dir.mkdir(parents=True, exist_ok=True)
        (self.template_dir / "clinical").mkdir(exist_ok=True)
        (self.template_dir / "research").mkdir(exist_ok=True)
        (self.template_dir / "patient").mkdir(exist_ok=True)
        
        # Initialize Jinja2
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=True
        )
    
    def generate_clinical_report(self, data: Dict[str, Any]) -> str:
        """
        Generate a clinical genomics report
        
        Args:
            data: Report data including patient info, variants, etc.
            
        Returns:
            HTML report as string
        """
        # Default clinical template
        template_html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Clinical Genomics Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; color: #333; }
        .header { border-bottom: 3px solid #2c3e50; padding-bottom: 20px; margin-bottom: 30px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 1px solid #ecf0f1; padding-bottom: 10px; }
        .patient-info { background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 20px 0; }
        .summary-box { background: #fff3cd; border: 1px solid #ffc107; padding: 15px; border-radius: 5px; }
        .pathogenic { background: #f8d7da; border: 1px solid #dc3545; padding: 10px; margin: 10px 0; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th { background: #34495e; color: white; padding: 10px; text-align: left; }
        td { padding: 10px; border-bottom: 1px solid #ecf0f1; }
        .footer { margin-top: 50px; padding-top: 20px; border-top: 1px solid #ecf0f1; color: #7f8c8d; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Clinical Genomics Report</h1>
        <p>Report Date: {{ report_date }}</p>
        <p>Report ID: {{ report_id }}</p>
    </div>

    <div class="patient-info">
        <h2>Patient Information</h2>
        <p><strong>Name:</strong> {{ patient_name }}</p>
        <p><strong>ID:</strong> {{ patient_id }}</p>
        <p><strong>Date of Birth:</strong> {{ patient_dob }}</p>
        <p><strong>Sample Type:</strong> {{ sample_type }}</p>
    </div>

    <div class="summary-box">
        <h2>Summary of Findings</h2>
        <p><strong>Total Variants Analyzed:</strong> {{ total_variants }}</p>
        <p><strong>Pathogenic Variants:</strong> {{ pathogenic_count }}</p>
        <p><strong>Likely Pathogenic:</strong> {{ likely_pathogenic_count }}</p>
        <p><strong>Variants of Uncertain Significance:</strong> {{ vus_count }}</p>
    </div>

    <h2>Clinically Significant Variants</h2>
    {% if pathogenic_variants %}
    <table>
        <thead>
            <tr>
                <th>Gene</th>
                <th>Variant</th>
                <th>Classification</th>
                <th>Condition</th>
                <th>Inheritance</th>
            </tr>
        </thead>
        <tbody>
            {% for variant in pathogenic_variants %}
            <tr>
                <td>{{ variant.gene }}</td>
                <td>{{ variant.change }}</td>
                <td>{{ variant.classification }}</td>
                <td>{{ variant.condition }}</td>
                <td>{{ variant.inheritance }}</td>
            </tr>
            {% endfor %}
        </tbody>
    </table>
    {% else %}
    <p>No pathogenic variants identified.</p>
    {% endif %}

    <h2>Recommendations</h2>
    <p>{{ recommendations }}</p>

    <div class="footer">
        <p>This report is for clinical use only and should be interpreted by qualified healthcare professionals.</p>
        <p>GenoScope Clinical Laboratory | {{ lab_certification }}</p>
    </div>
</body>
</html>"""
        
        # Create template from string
        from jinja2 import Template
        template = Template(template_html)
        
        # Prepare data with defaults
        report_data = {
            'report_date': datetime.now().strftime('%Y-%m-%d'),
            'report_id': data.get('report_id', 'GS-' + datetime.now().strftime('%Y%m%d%H%M%S')),
            'patient_name': data.get('patient', {}).get('name', 'Anonymous'),
            'patient_id': data.get('patient', {}).get('id', 'N/A'),
            'patient_dob': data.get('patient', {}).get('dob', 'N/A'),
            'sample_type': data.get('sample_type', 'Blood'),
            'total_variants': data.get('total_variants', 0),
            'pathogenic_count': data.get('pathogenic_count', 0),
            'likely_pathogenic_count': data.get('likely_pathogenic_count', 0),
            'vus_count': data.get('vus_count', 0),
            'pathogenic_variants': data.get('pathogenic_variants', []),
            'recommendations': data.get('recommendations', 'Genetic counseling recommended.'),
            'lab_certification': 'CLIA: 12D3456789 | CAP: 1234567'
        }
        
        return template.render(**report_data)
    
    def generate_research_report(self, data: Dict[str, Any]) -> str:
        """
        Generate a research genomics report with detailed statistics
        
        Args:
            data: Research data including variants, statistics, etc.
            
        Returns:
            HTML report as string
        """
        template_html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Research Genomics Report</title>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 30px; }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; }
        h1 { margin: 0; }
        .section { background: white; padding: 20px; margin: 20px 0; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; }
        .stat-card { background: #f8f9fa; padding: 15px; border-radius: 8px; }
        .stat-value { font-size: 2em; font-weight: bold; color: #667eea; }
        .stat-label { color: #6c757d; font-size: 0.9em; }
        table { width: 100%; border-collapse: collapse; }
        th { background: #f8f9fa; padding: 12px; text-align: left; }
        td { padding: 12px; border-bottom: 1px solid #dee2e6; }
        .chart-container { margin: 20px 0; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Research Genomics Analysis Report</h1>
        <p>Generated: {{ timestamp }}</p>
        <p>Analysis ID: {{ analysis_id }}</p>
    </div>

    <div class="section">
        <h2>Analysis Overview</h2>
        <div class="stat-grid">
            <div class="stat-card">
                <div class="stat-value">{{ total_variants }}</div>
                <div class="stat-label">Total Variants</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ snp_count }}</div>
                <div class="stat-label">SNPs</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ indel_count }}</div>
                <div class="stat-label">INDELs</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ mean_quality }}</div>
                <div class="stat-label">Mean Quality</div>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>Variant Distribution by Chromosome</h2>
        <table>
            <thead>
                <tr>
                    <th>Chromosome</th>
                    <th>Variant Count</th>
                    <th>Percentage</th>
                </tr>
            </thead>
            <tbody>
                {% for chr_data in chromosome_distribution %}
                <tr>
                    <td>{{ chr_data.chromosome }}</td>
                    <td>{{ chr_data.count }}</td>
                    <td>{{ chr_data.percentage }}%</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
    </div>

    <div class="section">
        <h2>Functional Impact Summary</h2>
        <ul>
            <li>Missense variants: {{ missense_count }}</li>
            <li>Nonsense variants: {{ nonsense_count }}</li>
            <li>Frameshift variants: {{ frameshift_count }}</li>
            <li>Splice site variants: {{ splice_count }}</li>
        </ul>
    </div>

    <div class="section">
        <h2>Population Frequency Analysis</h2>
        <p>Variants with MAF < 0.01: {{ rare_variants }}</p>
        <p>Variants with MAF 0.01-0.05: {{ low_freq_variants }}</p>
        <p>Common variants (MAF > 0.05): {{ common_variants }}</p>
    </div>

    <div class="section">
        <h2>Data Files</h2>
        <p>Input: {{ input_file }}</p>
        <p>Output: {{ output_files }}</p>
        <p>Log: {{ log_file }}</p>
    </div>
</body>
</html>"""
        
        from jinja2 import Template
        template = Template(template_html)
        
        # Prepare research data
        report_data = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'analysis_id': data.get('analysis_id', 'RES-' + datetime.now().strftime('%Y%m%d%H%M%S')),
            'total_variants': data.get('total_variants', 0),
            'snp_count': data.get('snp_count', 0),
            'indel_count': data.get('indel_count', 0),
            'mean_quality': data.get('mean_quality', 0),
            'chromosome_distribution': data.get('chromosome_distribution', []),
            'missense_count': data.get('missense_count', 0),
            'nonsense_count': data.get('nonsense_count', 0),
            'frameshift_count': data.get('frameshift_count', 0),
            'splice_count': data.get('splice_count', 0),
            'rare_variants': data.get('rare_variants', 0),
            'low_freq_variants': data.get('low_freq_variants', 0),
            'common_variants': data.get('common_variants', 0),
            'input_file': data.get('input_file', 'N/A'),
            'output_files': data.get('output_files', 'N/A'),
            'log_file': data.get('log_file', 'N/A')
        }
        
        return template.render(**report_data)
    
    def generate_patient_report(self, data: Dict[str, Any]) -> str:
        """
        Generate a patient-friendly genomics report
        
        Args:
            data: Simplified patient data
            
        Returns:
            HTML report as string
        """
        template_html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Your Genetic Test Results</title>
    <style>
        body { font-family: Georgia, serif; max-width: 800px; margin: 0 auto; padding: 40px 20px; line-height: 1.6; }
        .header { text-align: center; padding: 30px; background: #e8f5e9; border-radius: 10px; }
        h1 { color: #2e7d32; }
        .important { background: #fff3e0; padding: 20px; border-left: 4px solid #ff9800; margin: 20px 0; }
        .result-box { background: #f5f5f5; padding: 20px; margin: 20px 0; border-radius: 10px; }
        .positive { border-left: 4px solid #4caf50; }
        .negative { border-left: 4px solid #f44336; }
        .glossary { background: #e3f2fd; padding: 20px; border-radius: 10px; }
        .glossary dt { font-weight: bold; margin-top: 10px; }
        .contact { text-align: center; margin-top: 40px; padding: 20px; background: #f5f5f5; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Your Genetic Test Results</h1>
        <p>{{ patient_name }}</p>
        <p>Test Date: {{ test_date }}</p>
    </div>

    <div class="important">
        <h2>Important Information</h2>
        <p>This report contains your genetic test results. We recommend discussing these results with your healthcare provider or a genetic counselor who can help you understand what they mean for your health.</p>
    </div>

    <div class="result-box">
        <h2>Summary of Your Results</h2>
        <p>{{ summary }}</p>
        
        {% if key_findings %}
        <h3>Key Findings:</h3>
        <ul>
            {% for finding in key_findings %}
            <li>{{ finding }}</li>
            {% endfor %}
        </ul>
        {% endif %}
    </div>

    <div class="result-box">
        <h2>What This Means for You</h2>
        <p>{{ interpretation }}</p>
    </div>

    <div class="result-box">
        <h2>Next Steps</h2>
        <ol>
            {% for step in next_steps %}
            <li>{{ step }}</li>
            {% endfor %}
        </ol>
    </div>

    <div class="glossary">
        <h2>Understanding Genetic Terms</h2>
        <dl>
            <dt>Gene</dt>
            <dd>A segment of DNA that contains instructions for making proteins in your body.</dd>
            
            <dt>Variant</dt>
            <dd>A change in the DNA sequence compared to what is typically expected.</dd>
            
            <dt>Pathogenic</dt>
            <dd>A variant that is known to cause or increase risk for a genetic condition.</dd>
            
            <dt>Benign</dt>
            <dd>A variant that is not expected to cause health problems.</dd>
        </dl>
    </div>

    <div class="contact">
        <h3>Questions?</h3>
        <p>Contact your healthcare provider or our genetic counseling team:</p>
        <p>📞 1-800-GENETICS | 📧 counseling@genoscope.com</p>
    </div>
</body>
</html>"""
        
        from jinja2 import Template
        template = Template(template_html)
        
        # Prepare patient-friendly data
        report_data = {
            'patient_name': data.get('patient_name', 'Patient'),
            'test_date': data.get('test_date', datetime.now().strftime('%B %d, %Y')),
            'summary': data.get('summary', 'Your genetic test has been completed.'),
            'key_findings': data.get('key_findings', []),
            'interpretation': data.get('interpretation', 'Please discuss these results with your healthcare provider.'),
            'next_steps': data.get('next_steps', [
                'Schedule a follow-up appointment with your doctor',
                'Consider genetic counseling if recommended',
                'Share these results with family members if advised'
            ])
        }
        
        return template.render(**report_data)
    
    def save_report(self, html_content: str, output_path: str, format: str = 'html') -> str:
        """
        Save report to file
        
        Args:
            html_content: HTML report content
            output_path: Output file path
            format: Output format (html or pdf)
            
        Returns:
            Path to saved report
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if format == 'html':
            output_file.write_text(html_content)
            logger.info(f"Report saved to {output_file}")
            return str(output_file)
        
        elif format == 'pdf':
            # For PDF generation, we'd use weasyprint or similar
            # For now, just save as HTML with .pdf extension
            try:
                from weasyprint import HTML
                HTML(string=html_content).write_pdf(output_file)
                logger.info(f"PDF report saved to {output_file}")
            except ImportError:
                # Fallback to HTML if weasyprint not installed
                html_file = output_file.with_suffix('.html')
                html_file.write_text(html_content)
                logger.warning("PDF generation not available, saved as HTML")
                return str(html_file)
            
            return str(output_file)
        
        else:
            raise ValueError(f"Unsupported format: {format}")
