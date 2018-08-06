-- Alter DB to create a table that connects Transcripts to their GOs
-- Used to verify interest GOs appearance in the data set 

CREATE TABLE transcript_orf AS SELECT Transcript.transcript_id, ORF.orf_id FROM Transcript INNER JOIN ORF ON Transcript.transcript_id=ORF.transcript_id;

CREATE TABLE orf_blast AS SELECT ORF.orf_id, BlastDbase.FullAccession FROM ORF INNER JOIN BlastDbase ON ORF.orf_id=BlastDbase.TrinityID;

CREATE TABLE blast_uniprot AS SELECT BlastDbase.FullAccession, UniprotIndex.LinkID FROM BlastDbase INNER JOIN UniprotIndex ON BlastDbase.FullAccession=UniprotIndex.Accession;

CREATE TABLE uniprot_go AS SELECT go.id FROM UniprotIndex INNER JOIN go ON go.id=UniprotIndex.LinkID;

CREATE TABLE transcript_orf_blast AS SELECT transcript_orf.transcript_id, transcript_orf.orf_id, orf_blast.FullAccession FROM transcript_orf INNER JOIN orf_blast ON transcript_orf.orf_id=orf_blast.orf_id;

CREATE TABLE blast_uniprot_go AS SELECT DISTINCT blast_uniprot.FullAccession, uniprot_go.id FROM blast_uniprot INNER JOIN uniprot_go ON blast_uniprot.LinkID=uniprot_go.id;

-- Table final distinct contains relationship between transcripts and their GOs
CREATE TABLE final AS SELECT DISTINCT * FROM transcript_orf_blast INNER JOIN blast_uniprot_go ON transcript_orf_blast.FullAccession=blast_uniprot_go.FullAccession;

-- Remove middle (temporary) tables
DROP TABLE transcript_orf;
DROP TABLE orf_blast;
DROP TABLE blast_uniprot;
DROP TABLE uniprot_go;
DROP TABLE transcript_orf_blast;
DROP TABLE blast_uniprot_go;

