-- Alter DB to create a table that connects Transcripts to their GOs
-- Used to verify interest GOs appearance in the data set 
-- Used in SQLite

CREATE VIEW transcript_orf AS SELECT Transcript.transcript_id, ORF.orf_id FROM Transcript INNER JOIN ORF ON Transcript.transcript_id=ORF.transcript_id;

CREATE VIEW orf_blast AS SELECT ORF.orf_id, BlastDbase.FullAccession FROM ORF INNER JOIN BlastDbase ON ORF.orf_id=BlastDbase.TrinityID;

CREATE VIEW blast_uniprot AS SELECT BlastDbase.FullAccession, UniprotIndex.LinkID FROM BlastDbase INNER JOIN UniprotIndex ON BlastDbase.FullAccession=UniprotIndex.Accession;

CREATE VIEW uniprot_go AS SELECT go.id FROM UniprotIndex INNER JOIN go ON go.id=UniprotIndex.LinkID;

CREATE VIEW transcript_orf_blast AS SELECT transcript_orf.transcript_id, transcript_orf.orf_id, orf_blast.FullAccession FROM transcript_orf INNER JOIN orf_blast ON transcript_orf.orf_id=orf_blast.orf_id;

CREATE VIEW blast_uniprot_go AS SELECT DISTINCT blast_uniprot.FullAccession, uniprot_go.id FROM blast_uniprot INNER JOIN uniprot_go ON blast_uniprot.LinkID=uniprot_go.id;

-- Table final distinct contains relationship between transcripts and their GOs
CREATE VIEW final AS SELECT DISTINCT * FROM transcript_orf_blast INNER JOIN blast_uniprot_go ON transcript_orf_blast.FullAccession=blast_uniprot_go.FullAccession;

SELECT DISTINCT Transcript.sequence, ORF.peptide, BlastDbase.FullAccession, BlastDbase.GINumber, BlastDbase.PercentIdentity, BlastDbase.Evalue, UniprotIndex.Accession, UniprotIndex.AttributeType, go.id, go.name, go.namespace, go.def, TaxonomyIndex.NCBITaxonomyAccession, TaxonomyIndex.TaxonomyValue, PFAMreference.pfam_accession, PFAMreference.pfam_domainname, PFAMreference.pfam_domaindescription
FROM Transcript 
INNER JOIN ORF ON ORF.transcript_id=Transcript.transcript_id
INNER JOIN BlastDbase ON (BlastDbase.TrinityID=ORF.orf_id OR BlastDbase.TrinityID=ORF.transcript_id)
INNER JOIN UniprotIndex ON  BlastDbase.FullAccession=UniprotIndex.Accession
INNER JOIN go ON go.id=UniprotIndex.LinkID
LEFT JOIN TaxonomyIndex ON TaxonomyIndex.NCBITaxonomyAccession=UniprotIndex.LinkId
INNER JOIN pfam2go ON go.id=pfam2go.go_id 
LEFT JOIN PFAMreference ON PFAMreference.pfam_accession=pfam2go.pfam_acc;

-- clean slate
DROP VIEW transcript_orf;
DROP VIEW orf_blast;
DROP VIEW blast_uniprot;
DROP VIEW uniprot_go;
DROP VIEW transcript_orf_blast;
DROP VIEW blast_uniprot_go;
DROP VIEW final;

-- core data
CREATE VIEW core_data AS SELECT DISTINCT Transcript.sequence, ORF.peptide, BlastDbase.FullAccession, BlastDbase.GINumber, BlastDbase.PercentIdentity, BlastDbase.Evalue, UniprotIndex.Accession, UniprotIndex.AttributeType, go.id, go.name, go.namespace, go.def, UniprotIndex.LinkID, HMMERDbase.HMMERDomain, HMMERDbase.HMMERTDomainDescription, HMMERDbase.FullSeqScore, PFAMreference.pfam_accession, PFAMreference.pfam_domainname, PFAMreference.pfam_domaindescription, TaxonomyIndex.NCBITaxonomyAccession, TaxonomyIndex.TaxonomyValue, eggNOGIndex.eggNOGIndexTerm, eggNOGIndex.eggNOGDescriptionValue
FROM Transcript 
INNER JOIN ORF ON ORF.transcript_id=Transcript.transcript_id
INNER JOIN BlastDbase ON (BlastDbase.TrinityID=ORF.orf_id OR BlastDbase.TrinityID=ORF.transcript_id)
INNER JOIN UniprotIndex ON  BlastDbase.FullAccession=UniprotIndex.Accession
INNER JOIN go ON go.id=UniprotIndex.LinkID
INNER JOIN HMMERDbase ON (HMMERDbase.QueryProtID=ORF.orf_id OR HMMERDbase.QueryProtID=ORF.transcript_id)
INNER JOIN PFAMreference ON HMMERDbase.pfam_id=PFAMreference.pfam_accession
LEFT JOIN TaxonomyIndex ON UniprotIndex.LinkId=TaxonomyIndex.NCBITaxonomyAccession
LEFT JOIN eggNOGIndex ON eggNOGIndex.eggNOGIndexTerm=UniprotIndex.LinkId;


