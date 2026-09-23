# Raw-count doublet investigation

## Provenance
- DATASET-DERIVED: all metrics below use `ad.layers['counts']`, the supplied raw UMI count layer.
- TRANSCRIPT-DERIVED: programme genes are the owner-supplied clues; PPBP/PF4 are explicitly platelet-contamination clues.
- MODEL PRETRAINED KNOWLEDGE: high depth, incompatible programme co-expression, and doublet-style interpretation are general biological/technical knowledge, not statements from the transcript.

## Method
- Calculated raw total counts and raw detected genes per cell.
- Calculated robust z-scores using the dataset-wide median and MAD.
- For each supplied programme, counted raw positive marker genes and raw counts.
- Flagged a possible mixed cell when at least two programme groups had multiple positive markers; the platelet clue required one positive gene. This is a screening rule, not a validated doublet classifier.

## Cluster summary
cluster  n_cells  median_raw_counts  max_raw_counts  median_raw_genes  max_raw_genes  mixed_cells  median_programmes  max_programmes
      0     1197             2310.0          6329.0             808.0           1815          260                1.0               4
      1      489             2305.0          5682.0             848.0           1649           93                1.0               4
      2      445             1954.0          5343.0             826.0           1868          143                1.0               3
      3      347             1759.0          8011.0             671.0           1996           61                1.0               3
      4      163             3781.0          7167.0            1262.0           1934           39                1.0               3
      5       36             5293.0          8415.0            1565.5           2013           13                1.0               3
      6       13              916.0          8875.0             345.0           2413            4                1.0               3
      7       10             8498.5         15818.0            2354.0           3400            8                2.0               3

## Focus cells
                 cluster  raw_total_counts  raw_n_genes  raw_counts_robust_z  raw_genes_robust_z  n_strong_programmes  possible_mixed_programmes
index                                                                                                                                           
GTGATGACAAGTGA-1       0            6329.0         1815             5.575341            5.143636                    3                       True
ACTCTCCTGCATAC-1       0            5850.0         1739             4.929179            4.752328                    0                      False
TCGGACCTGTACAC-1       0            5599.0         1566             4.590584            3.861588                    2                       True
CTAAGGACCGTTAG-1       0            5484.0         1425             4.435452            3.135610                    1                      False
GGCATATGTGTGAC-1       0            5303.0         1584             4.191286            3.954266                    2                       True
TACGTACTACGGAG-1       0            5280.0         1606             4.160259            4.067540                    4                       True
TTCAGTACTCAAGC-1       0            5218.0         1377             4.076622            2.888468                    1                      False
AGCTGCCTTTCATC-1       0            5212.0         1703             4.068528            4.566972                    2                       True
TGACGCCTGTACCA-1       0            5155.0         1646             3.991636            4.273491                    2                       True
CGATCAGAGGTACT-1       0            5139.0         1338             3.970053            2.687665                    2                       True
TGAGACACAAGGTA-1       0            5135.0         1546             3.964657            3.758613                    1                      False
ATGATATGAAACAG-1       0            4921.0         1399             3.675975            3.001741                    3                       True
ATGTAAACTCTCCG-1       0            4766.0         1323             3.466882            2.610434                    2                       True
CTATACTGTTCGTT-1       0            4760.0         1520             3.458789            3.624744                    1                      False
GAACCTGAACGTGT-1       0            4754.0         1489             3.450695            3.465132                    2                       True
GATCGAACCGAGAG-1       0            4587.0         1425             3.225415            3.135610                    2                       True
GTATTAGAAACAGA-1       0            3487.0         1426             1.741535            3.140758                    1                      False
AACCTACTGTGAGG-1       1            5682.0         1649             4.702549            4.288937                    1                      False
ACCCAGCTGTTAGC-1       1            5534.0         1546             4.502901            3.758613                    2                       True
CTAGAGACTTTGGG-1       1            5515.0         1621             4.477270            4.144771                    2                       True
GAACAGCTAACTGC-1       1            5489.0         1631             4.442196            4.196259                    1                      False
AGAGCGGAGGCAAG-1       1            5226.0         1422             4.087414            3.120163                    1                      False
TAAGCGTGGACAAA-1       1            5189.0         1634             4.037502            4.211706                    2                       True
CCGAAAACCTTGTT-1       1            5144.0         1512             3.976798            3.583554                    1                      False
GTTAGGTGCCAGTA-1       1            4971.0         1496             3.743424            3.501173                    2                       True
ATGCACGATTGGTG-1       1            4854.0         1551             3.585593            3.784357                    1                      False
CCCAGACTGCCTTC-1       1            4754.0         1378             3.450695            2.893617                    1                      False
GTCCCATGTGGTGT-1       1            4744.0         1467             3.437205            3.351859                    2                       True
GGTCTAGATAGCGT-1       1            4660.0         1460             3.323890            3.315817                    1                      False
AAAGCAGATATCGG-1       1            4584.0         1422             3.221368            3.120163                    1                      False
GTAACGTGACCTCC-1       1            4521.0         1321             3.136382            2.600136                    1                      False
ATACTCTGCTTCGC-1       1            4471.0         1331             3.068933            2.651624                    2                       True
CCTCTACTCTTCGC-1       1            4464.0         1322             3.059490            2.605285                    2                       True
CAGTTTACACACGT-1       2            5343.0         1868             4.245245            5.416521                    2                       True
GTGTATCTAGTAGA-1       2            4442.0         1426             3.029813            3.140758                    2                       True
GGCGACACTGCCCT-1       2            4438.0         1494             3.024417            3.490876                    1                      False
TAATGCCTCGTCTC-1       2            4211.0         1667             2.718198            4.381616                    2                       True
TATGGGTGCATCAG-1       2            3544.0         1457             1.818427            3.300371                    2                       True
GCGCGATGGTGCAT-1       2            3225.0         1410             1.388102            3.058378                    1                      False
CAGGTTGAGGATCT-1       3            8011.0         1996             7.844327            6.075566                    2                       True
ACGAAGCTCTGAGT-1       3            7062.0         1852             6.564144            5.334141                    1                      False
AAGCACTGGTTCTT-1       3            6153.0         1713             5.337920            4.618460                    2                       True
CTAGGATGAGCCTA-1       3            5867.0         1741             4.952111            4.762626                    1                      False
TACTGTTGCTGAAC-1       3            5187.0         1650             4.034804            4.294086                    2                       True
AAAGGCCTGTCTAG-1       3            4973.0         1445             3.746122            3.238585                    3                       True
AAACATTGAGCTAC-1       3            4903.0         1352             3.651693            2.759748                    1                      False
CTATGTACTGTTTC-1       3            4824.0         1481             3.545124            3.423942                    2                       True
CATACTTGGGTTAC-1       4            7167.0         1934             6.705787            5.756341                    1                      False
GAGCATACTTTGCT-1       4            6691.0         1750             6.063672            4.808965                    2                       True
ATTTAGGAACCATG-1       4            5677.0         1662             4.695805            4.355872                    1                      False
AAATCAACCCTATT-1       4            5676.0         1541             4.694456            3.732869                    3                       True
GCTCAAGAACCATG-1       4            5552.0         1674             4.527182            4.417657                    2                       True
TGCGTAGATGGTCA-1       4            5486.0         1557             4.438149            3.815249                    1                      False
CATGCGCTAGTCAC-1       4            5344.0         1659             4.246594            4.340425                    1                      False
GAAACCTGGACTAC-1       4            5326.0         1690             4.222312            4.500038                    3                       True
TACGCCACTCCCAC-1       4            5286.0         1618             4.168353            4.129325                    2                       True
CATCAGGATCCTAT-1       4            5281.0         1582             4.161608            3.943969                    1                      False
AGTCTACTTGCATG-1       4            5201.0         1569             4.053689            3.877035                    1                      False
ATCTGTTGCCTTCG-1       4            5197.0         1550             4.048294            3.779208                    1                      False
GACAGTACTTCGGA-1       4            5193.0         1509             4.042898            3.568108                    1                      False
CGCAGGACTTGTCT-1       4            5162.0         1687             4.001079            4.484591                    2                       True
GTAGCTGAAGCTAC-1       4            5102.0         1526             3.920140            3.655637                    2                       True
GTACCCTGACAGTC-1       4            5049.0         1548             3.848644            3.768910                    1                      False
TAAGCGTGTGCTCC-1       4            4972.0         1502             3.744773            3.532066                    1                      False
GACAGTACGAGCTT-1       4            4929.0         1546             3.686767            3.758613                    2                       True
TATAAGTGTGGTGT-1       4            4913.0         1517             3.665183            3.609298                    1                      False
TGTACTTGCTCTAT-1       4            4907.0         1479             3.657089            3.413644                    1                      False
ATTGCACTTGCTTT-1       4            4888.0         1519             3.631458            3.619595                    2                       True
CTTCACCTACCTGA-1       4            4888.0         1528             3.631458            3.665935                    2                       True
AGCCTCACTGTCAG-1       4            4811.0         1393             3.527587            2.970849                    1                      False
CATGCGCTTTGCAG-1       4            4790.0         1515             3.499258            3.599000                    1                      False
CCACCATGGACGAG-1       4            4779.0         1558             3.484419            3.820398                    1                      False
TGAGGTACGAACCT-1       4            4752.0         1521             3.447997            3.629893                    1                      False
AACATTGATGGGAG-1       4            4711.0         1458             3.392689            3.305520                    2                       True
AGTGACTGCAACTG-1       4            4696.0         1412             3.372454            3.068676                    1                      False
GAAAGATGCTGATG-1       4            4688.0         1464             3.361662            3.336412                    1                      False
TCATTCGATACAGC-1       4            4629.0         1416             3.282072            3.089271                    1                      False
AACCTTTGTACGCA-1       4            4617.0         1461             3.265884            3.320966                    2                       True
TACATCACTGAACC-1       4            4600.0         1507             3.242952            3.557810                    1                      False
TCACCGTGCTCGCT-1       4            4572.0         1477             3.205180            3.403347                    2                       True
GTGGATTGCACTAG-1       4            4560.0         1466             3.188992            3.346710                    1                      False
CAACGATGCGCAAT-1       4            4517.0         1407             3.130986            3.042932                    1                      False
CTCAGGCTCGTTGA-1       4            4502.0         1463             3.110751            3.331264                    1                      False
AGAGTCTGGTCGTA-1       4            4469.0         1496             3.066235            3.501173                    1                      False
AGTATAACTTGTCT-1       4            4454.0         1498             3.046000            3.511471                    1                      False
GTTGGATGTTTACC-1       4            4436.0         1362             3.021719            2.811236                    1                      False
ACAGCAACCTCAAG-1       4            4431.0         1472             3.014974            3.377603                    1                      False
GATATTGACGAGTT-1       4            4409.0         1446             2.985296            3.243734                    2                       True
ATCGCAGAATCTCT-1       4            4337.0         1418             2.888170            3.099568                    2                       True
GGGAACGAAGCTCA-1       4            4319.0         1554             2.863888            3.799803                    1                      False
TGATCGGACTGACA-1       4            4211.0         1433             2.718198            3.176800                    2                       True
GAGTTGTGTATGCG-1       4            4175.0         1416             2.669635            3.089271                    1                      False
CTGAAGTGAAGCCT-1       4            4070.0         1434             2.527991            3.181949                    1                      False
GAAGTAGATCCAAG-1       4            4035.0         1408             2.480777            3.048080                    1                      False
GGGCCAACCTTGGA-1       5            8415.0         2013             8.389317            6.163095                    0                      False
ACGAGGGACAGGAG-1       5            7928.0         1991             7.732362            6.049822                    1                      False
AAGCCATGAACTGC-1       5            7064.0         1871             6.566842            5.431968                    1                      False
CGATCAGATGTGAC-1       5            6908.0         1974             6.356401            5.962292                    1                      False
ATACCGGAATGCTG-1       5            6580.0         1859             5.913935            5.370182                    1                      False
TTGAGGACTACGCA-1       5            6342.0         1792             5.592877            5.025214                    1                      False
TGCAATCTTCAGGT-1       5            6316.0         1902             5.557804            5.591580                    1                      False
TGTAGGTGCTCTAT-1       5            6280.0         1906             5.509241            5.612175                    2                       True
GACATTCTCCACCT-1       5            6191.0         1749             5.389181            4.803816                    1                      False
GAAAGTGAAAGTGA-1       5            6175.0         1571             5.367598            3.887332                    1                      False
TTATGGCTTATGGC-1       5            6164.0         1782             5.352759            4.973726                    1                      False
GTTAACCTTGCTTT-1       5            6097.0         1643             5.262377            4.258045                    2                       True
ACCCGTTGCTTCTA-1       5            6083.0         1852             5.243491            5.334141                    3                       True
GCGTAAACACGGTT-1       5            5924.0         1625             5.029003            4.165367                    3                       True
TTTAGCTGTACTCT-1       5            5671.0         1560             4.687711            3.830696                    2                       True
GAAGGTCTTAAAGG-1       5            5637.0         1627             4.641846            4.175664                    1                      False
ATTGTAGATTCCCG-1       5            5496.0         1647             4.451639            4.278640                    1                      False
ACGTGATGCCATGA-1       5            5437.0         1414             4.372049            3.078973                    1                      False
AGCACTGATGCTTT-1       5            5149.0         1605             3.983542            4.062391                    2                       True
TTACTCGACGCAAT-1       5            5030.0         1602             3.823014            4.046945                    1                      False
CTCCACGAGAGATA-1       5            4633.0         1445             3.287468            3.238585                    2                       True
CGCCTAACGAATGA-1       5            4597.0         1455             3.238905            3.290073                    1                      False
AATTACGAATTCCT-1       5            4510.0         1313             3.121543            2.558946                    1                      False
GGACCGTGGGAACG-1       5            4455.0         1425             3.047349            3.135610                    1                      False
TTTCGAACACCTGA-1       5            4455.0         1539             3.047349            3.722571                    0                      False
AATGCGTGGACGGA-1       5            4432.0         1266             3.016323            2.316953                    2                       True
ACGAACTGGCTATG-1       6            8875.0         2413             9.009848            8.222609                    3                       True
GCGCATCTGGTTAC-1       6            2513.0          979             0.427627            0.839252                    2                       True
GAGTTGTGGTAGCT-1       6            2405.0          849             0.281937            0.169910                    1                      False
TAACACCTTGTTTC-1       6             999.0          394            -1.614731           -2.172787                    1                      False
GTCATACTTCGCCT-1       6             979.0          359            -1.641711           -2.352994                    1                      False
TTACGTACGTTCAG-1       6             942.0          321            -1.691623           -2.548648                    1                      False
ATTCAGCTCATTGG-1       6             916.0          345            -1.726696           -2.425077                    2                       True
ATCATCTGACACCA-1       6             736.0          366            -1.969513           -2.316953                    1                      False
GACGCTCTCTCTCG-1       6             682.0          266            -2.042358           -2.831831                    1                      False
AGTCTTACTTCGGA-1       6             663.0          273            -2.067989           -2.795790                    1                      False
GGAACACTTCAGAC-1       6             623.0          278            -2.121948           -2.770046                    1                      False
GGCATATGGGGAGT-1       6             575.0          212            -2.186699           -3.109866                    2                       True
ACCCACTGGTTCAG-1       6             556.0          338            -2.212330           -2.461119                    1                      False
CCAGTCTGCGGAGA-1       7           15818.0         3400            18.375826           13.304459                    3                       True
TTACTCGAACGTTG-1       7           15287.0         3302            17.659517           12.799878                    3                       True
AGAGGTCTACAGCT-1       7           10762.0         2798            11.555376           10.204891                    2                       True
GGCACGTGTGAGAA-1       7           10349.0         2685            10.998246            9.623078                    3                       True
GCGAAGGAGAGCTT-1       7           10275.0         2703            10.898422            9.715756                    2                       True
CGATACGACAGGAG-1       7            6722.0         2023             6.105491            6.214583                    1                      False
GCCTCAACTCTTTG-1       7            5981.0         1954             5.105895            5.859317                    2                       True
ACTTCAACAAGCAA-1       7            5416.0         1777             4.343720            4.947982                    1                      False
CACCGGGACTTCTA-1       7            4295.0         1610             2.831512            4.088135                    2                       True
ACGTCGCTCCTGAA-1       7            3924.0         1589             2.331040            3.980010                    2                       True

## Assessment
- Cluster 6: the cluster contains 13 cells and includes a high-depth cell; its transcript-supported PP4/PF4 platelet-contamination signal should be investigated. The raw-count screen alone does not establish a doublet.
- Cluster 7: the cluster has high depth and a few cells may show multiple programme signals, but broad marker positivity is not by itself proof of a doublet. Its suspected cycling interpretation is model knowledge, not transcript knowledge.
- A stronger doublet call would require per-cell incompatible raw-count co-expression plus supporting evidence from a dedicated method such as Scrublet and review of the flagged cell’s neighbourhood.

## Outputs
- `results/doublet_raw_cell_metrics.csv`
- `results/doublet_cluster_summary.csv`
- `results/doublet_focus_cells.csv`
