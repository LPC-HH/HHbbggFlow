#!/bin/bash

# Define the base directory
BASE_DIR="/eos/uscms/store/group/lpcdihiggsboost/tsievert/HiggsDNA_parquet/v1"

# Define output JSON file
OUTPUT_FILE="output.json"

# Declare associative arrays for luminosity, cross-sections, and branching fractions
declare -A LUMI_VALUES=(
    ["Run3_2022preEE"]="7.9804",
    ["Run3_2022postEE"]="26.6717"
    # Add more ERA values with their luminosities here
)

declare -A XS_VALUES=(
    #in fb
    ["ttHToGG"]=506.5,
    ["GluGluToHH"]=34.43,
    ["GGJets"]=88750,
    ["GJetPt20To40"]=242500,
    ["GJetPt40"]=919100,
    ["GluGluHToGG"]=48520,
    ["VBFHToGG"]=3779,
    ["VHToGG"]=2251.4, #1369 + 882.4
    ["Data_EraC"]=1,
    ["Data_EraD"]=1,
    ["Data_EraE"]=1,
    ["Data_EraF"]=1,
    ["Data_EraG"]=1
    # Add more SAMPLE values with their cross-section values here
)

declare -A BF_VALUES=(
    ["ttHToGG"]=0.00228, #htogg
    ["GluGluToHH"]=0.00265, #hh to bbgg
    ["GGJets"]=1,
    ["GJetPt20To40"]=1,
    ["GJetPt40"]=1,
    ["GluGluHToGG"]=0.00228,
    ["VBFHToGG"]=0.00228,
    ["VHToGG"]=0.00228,
    ["Data_EraC"]=1,
    ["Data_EraD"]=1,
    ["Data_EraE"]=1,
    ["Data_EraF"]=1,
    ["Data_EraG"]=1
    # Add more SAMPLE values with their branching fraction values here
)

# Start building the JSON structure
echo "{" > $OUTPUT_FILE
echo '    "data_eras" : {' >> $OUTPUT_FILE

# Loop through each ERA directory
for ERA in "${!LUMI_VALUES[@]}"; do
    ERA_DIR="$BASE_DIR/$ERA"
    
    if [ -d "$ERA_DIR" ]; then
        echo '        "'$ERA'" : {' >> $OUTPUT_FILE
        echo '            "lumi" : "'${LUMI_VALUES[$ERA]}'",' >> $OUTPUT_FILE
        echo '            "file_prefix" : "'$BASE_DIR'",' >> $OUTPUT_FILE
        echo '            "samples" : {' >> $OUTPUT_FILE

        # Loop through each SAMPLE directory under ERA
        for SAMPLE in "${!XS_VALUES[@]}"; do
            SAMPLE_DIR="$ERA_DIR/$SAMPLE/nominal"

            if [ -d "$SAMPLE_DIR" ]; then
                # Get the list of files
                FILES=$(ls "$SAMPLE_DIR")

                # Start the sample's JSON entry
                echo '                "'$SAMPLE'" : {' >> $OUTPUT_FILE
                echo '                    "xs" : '${XS_VALUES[$SAMPLE]}',' >> $OUTPUT_FILE
                echo '                    "bf" : '${BF_VALUES[$SAMPLE]}',' >> $OUTPUT_FILE
                echo '                    "files" : {' >> $OUTPUT_FILE
                echo '                        "nominal" : [' >> $OUTPUT_FILE

                # Loop through each file and add it to the JSON structure
                for FILE in $FILES; do
                    echo '                            "'$FILE'",' >> $OUTPUT_FILE
                done

                # Remove the trailing comma from the last file entry
                sed -i '$ s/,$//' $OUTPUT_FILE

                # Close the sample's JSON entry
                echo '                        ]' >> $OUTPUT_FILE
                echo '                    }' >> $OUTPUT_FILE
                echo '                },' >> $OUTPUT_FILE
            fi
        done

        # Remove the trailing comma from the last sample entry
        sed -i '$ s/,$//' $OUTPUT_FILE

        # Close the ERA's JSON entry
        echo '            }' >> $OUTPUT_FILE
        echo '        },' >> $OUTPUT_FILE
    fi
done

# Remove the trailing comma from the last ERA entry
sed -i '$ s/,$//' $OUTPUT_FILE

# Close the JSON structure
echo '    }' >> $OUTPUT_FILE
echo "}" >> $OUTPUT_FILE

# Notify the user
echo "JSON output written to $OUTPUT_FILE"
