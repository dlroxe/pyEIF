This directory has various data files for evaluating Entrez-related logic.

1. first_20_Hs.data has the first 20 records from Hs.data.

2. first_3_short_Hs.data has the first 3 records form Hs.Data, stripped of all
   but the most essential rows.

3. short_Hs.data has a record for EVERY record in Hs.data, stripped of everything
   except ID, GENE, and GENE_ID, along with 'SCOUNT      0' (to force the
   BioPython parser to accept that there are no 'SEQUENCE' lines).
