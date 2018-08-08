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

