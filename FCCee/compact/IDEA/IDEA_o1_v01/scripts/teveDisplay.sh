### Script to look at detector geometry (warning: very slow for complex detectors)
### and over bad connection

xml=$1

if [[ -n "$xml" ]]; then # test to see if not empty
    teveDisplay -ui -compact ${xml}
else
    echo "argument error, please provide an xml file as input argument!"
fi

