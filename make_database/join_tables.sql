-- Alter DB to create a table that connects Transcripts to their GOs
-- Used to verify interest GOs appearance in the data set 
-- Used in SQLite

SELECT 
    Transcript.sequence AS Sequence, 
    ORF.peptide AS Peptide, 
    Transcript.origin_pair AS `From Pair`,
    BlastDbase.FullAccession AS `Blast Acession`, 
    BlastDbase.PercentIdentity AS `Blast Identity`, 
    BlastDbase.Evalue AS `Blast E-value`, 
    UniprotIndex.Accession AS `Uniprot Acession`, 
    UniprotIndex.AttributeType AS `Uniprot Attr Type`, 
    go.id AS `GO`, 
    go.name AS `GO Name`, 
    go.namespace AS `GO Namespace`, 
    go.def AS `GO Def`, 
    TaxonomyIndex.NCBITaxonomyAccession AS `Taxonomy Index`, 
    TaxonomyIndex.TaxonomyValue AS `Taxonomy Value`, 
    PFAMreference.pfam_accession, 
    PFAMreference.pfam_domainname, 
    PFAMreference.pfam_domaindescription
FROM 
    Transcript 
LEFT JOIN ORF ON 
    ORF.transcript_id=Transcript.transcript_id
LEFT JOIN BlastDbase ON 
    (BlastDbase.TrinityID=ORF.orf_id OR BlastDbase.TrinityID=Transcript.transcript_id)
LEFT JOIN UniprotIndex ON 
    BlastDbase.FullAccession=UniprotIndex.Accession
LEFT JOIN go ON 
    go.id=UniprotIndex.LinkID
LEFT JOIN TaxonomyIndex ON 
    TaxonomyIndex.NCBITaxonomyAccession=UniprotIndex.LinkId
LEFT JOIN pfam2go ON 
    go.id=pfam2go.go_id 
LEFT JOIN PFAMreference ON 
    PFAMreference.pfam_accession=pfam2go.pfam_acc;


