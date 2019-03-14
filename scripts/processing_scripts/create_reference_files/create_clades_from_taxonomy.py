import pandas
import sys

#### Example input file (tab-delimited).
# assembly        Superkingdom    Phylum  Class   Order   Family  Genus   Species
# 637000057       Bacteria        Proteobacteria  Gammaproteobacteria     Enterobacterales        Enterobacteriaceae      Candidatus Blochmannia  Candidatus Blochmannia pennsylvanicus
# 637000058       Bacteria        Proteobacteria  Alphaproteobacteria     Pelagibacterales        Pelagibacteraceae       Candidatus Pelagibacter Candidatus Pelagibacter ubique
# 637000059       Bacteria        Chlamydiae      Chlamydiia      Parachlamydiales        Parachlamydiaceae       Candidatus Protochlamydia       Candidatus Protochlamydia amoebophila


if len(sys.argv) != 2:
    sys.exit("Need taxonomy table as input. See top of file for example of what this looks like.")

taxa = pandas.read_csv(sys.argv[1], sep="\t")

# Drop rows with NA values.
taxa.dropna(inplace=True)

# Combine all taxa columns at each level.
sep = "|"

taxa['Phylum_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum']

taxa['Class_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum'] + sep + \
                        taxa['Class']

taxa['Order_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum'] + sep + \
                        taxa['Class'] + sep + taxa['Order']

taxa['Family_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum'] + sep + \
                        taxa['Class'] + sep + taxa['Order'] + sep + \
                        taxa['Family']

taxa['Genus_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum'] + sep + \
                        taxa['Class'] + sep + taxa['Order'] + sep + \
                        taxa['Family'] + sep + taxa['Genus']

taxa['Species_lineage'] = taxa['Superkingdom'] + sep + taxa['Phylum'] + sep + \
                        taxa['Class'] + sep + taxa['Order'] + sep + \
                        taxa['Family'] + sep + taxa['Genus'] +  sep + \
                        taxa['Species']

lineages = ['Species_lineage', 'Genus_lineage', 'Family_lineage',
            'Order_lineage', 'Class_lineage', 'Phylum_lineage', 'Superkingdom']

lineage_groupings = {}

past_taxa = set()

for i in range(len(lineages)):

    current_lineage = lineages[i]

    if i == 0:
        past_lineage = "assembly"
    else:
        past_lineage  = lineages[i - 1]

    # Get all unique values for this lineage column.
    unique_taxa = list(taxa[current_lineage].unique())

    for taxon in unique_taxa:
        taxon_descendents = set(taxa.loc[taxa[current_lineage] == taxon, past_lineage])

        if len(taxon_descendents.intersection(past_taxa)) > 0:
            sys.exit("Error - descendents already seen.")

        # Make a new set with these added descendents.
        past_taxa = past_taxa.union(taxon_descendents)

        # Need to get clade format unless dealing with assembly ids.
        if i > 0:
            clade_descendents = [lineage_groupings[descendent] for descendent in taxon_descendents]
            taxon_descendents = list(clade_descendents)
        else:
            taxon_descendents = [str(descendent) for descendent in taxon_descendents]

        # Write descendents in the newick tree clade format representing that
        # taxon.
        if len(taxon_descendents) > 1:
            clade_format = "(" + ", ".join(taxon_descendents) + ")"
        else:
            clade_format = taxon_descendents[0]

        lineage_groupings[taxon] = clade_format


# Print out newick format.
unique_Superkingdom = list(taxa["Superkingdom"].unique())

Superkingdom_clades = [lineage_groupings[kingdom] for kingdom in unique_Superkingdom]

if len(unique_Superkingdom) == 1:
    print(Superkingdom_clades[0] + ";")
else:
    print("(" + ", ".join(Superkingdom_clades) + ");")
