-- Contact-based interactions
SELECT e1.ensembl_exon_id as exon1,
    e2.ensembl_exon_id as exon2,
    p1.uniprot_id as protein1,
    p2.uniprot_id as protein2,
    ei.pdb_id,
    ei.jaccard_percent,
    ei.exon1_coverage_percent,
    ei.exon2_coverage_percent
FROM eei_interactions ei
    JOIN exons e1 ON ei.exon1_id = e1.exon_id
    JOIN exons e2 ON ei.exon2_id = e2.exon_id
    JOIN proteins p1 ON ei.protein1_id = p1.protein_id
    JOIN proteins p2 ON ei.protein2_id = p2.protein_id
    JOIN eei_methods em ON ei.method_id = em.method_id
WHERE em.method_name = 'contact_based'
ORDER BY ei.jaccard_percent DESC
LIMIT 10;
SELECT *
FROM get_eei_statistics();
-- PISA interactions
SELECT COUNT(*) as total_pisa_eeis,
    COUNT(DISTINCT ei.exon1_id) as unique_exon1,
    COUNT(DISTINCT ei.exon2_id) as unique_exon2,
    COUNT(DISTINCT ei.protein1_id) as unique_proteins1,
    COUNT(DISTINCT ei.protein2_id) as unique_proteins2,
    COUNT(DISTINCT ei.pdb_id) as unique_pdb_structures
FROM eei_interactions ei
    JOIN eei_methods em ON ei.method_id = em.method_id
WHERE em.method_name = 'PISA';
-- Overall database statistics
SELECT *
FROM get_eei_statistics();
-- Compare contact-based vs PISA
SELECT em.method_name,
    COUNT(*) as total_eeis,
    COUNT(DISTINCT ei.exon1_id) as unique_exon1,
    COUNT(DISTINCT ei.exon2_id) as unique_exon2,
    COUNT(DISTINCT ei.pdb_id) as unique_pdb_structures
FROM eei_interactions ei
    JOIN eei_methods em ON ei.method_id = em.method_id
WHERE em.method_name IN ('contact_based', 'PISA')
GROUP BY em.method_name
ORDER BY em.method_name;
-- Sample some PISA interactions with best (most negative) free energy
SELECT e1.ensembl_exon_id as exon1,
    e2.ensembl_exon_id as exon2,
    p1.uniprot_id as protein1,
    p2.uniprot_id as protein2,
    ei.pdb_id,
    ROUND(epa.free_energy::numeric, 2) as free_energy,
    epa.hydrogen_bonds,
    epa.salt_bridges
FROM eei_interactions ei
    JOIN exons e1 ON ei.exon1_id = e1.exon_id
    JOIN exons e2 ON ei.exon2_id = e2.exon_id
    JOIN proteins p1 ON ei.protein1_id = p1.protein_id
    JOIN proteins p2 ON ei.protein2_id = p2.protein_id
    JOIN eei_methods em ON ei.method_id = em.method_id
    JOIN eei_pisa_attributes epa ON ei.eei_id = epa.eei_id
WHERE em.method_name = 'PISA'
ORDER BY epa.free_energy ASC -- Most stable interactions first
LIMIT 10;
-- EPPIC
SELECT COUNT(*) as total_eppic_eeis,
    COUNT(DISTINCT ei.exon1_id) as unique_exon1,
    COUNT(DISTINCT ei.exon2_id) as unique_exon2,
    COUNT(DISTINCT ei.protein1_id) as unique_proteins1,
    COUNT(DISTINCT ei.protein2_id) as unique_proteins2,
    COUNT(DISTINCT ei.pdb_id) as unique_pdb_structures
FROM eei_interactions ei
    JOIN eei_methods em ON ei.method_id = em.method_id
WHERE em.method_name = 'EPPIC';
SELECT COUNT(*) as total_with_orthology,
    AVG(eom.confidence) as avg_confidence,
    MIN(eom.confidence) as min_confidence,
    MAX(eom.confidence) as max_confidence,
    AVG(eom.identity1) as avg_identity1,
    AVG(eom.identity2) as avg_identity2
FROM eei_orthology_mapping eom
    JOIN eei_interactions ei ON eom.eei_id = ei.eei_id
    JOIN eei_methods em ON ei.method_id = em.method_id
WHERE em.method_name = 'predicted_contact';
-- Compare all methods (experimental + predicted)
SELECT em.method_name,
    em.method_type,
    COUNT(*) as total_eeis,
    COUNT(DISTINCT ei.exon1_id) as unique_exon1,
    COUNT(DISTINCT ei.exon2_id) as unique_exon2,
    COUNT(DISTINCT ei.pdb_id) as unique_pdb_structures
FROM eei_interactions ei
    JOIN eei_methods em ON ei.method_id = em.method_id
GROUP BY em.method_name,
    em.method_type
ORDER BY em.method_type,
    em.method_name;