#!/bin/bash

WORK_DIR=$(pwd)
SIRIL_PATH=/Applications/Siril.app/Contents/MacOS/siril-cli

#--------------------------------------------
# Functions
#--------------------------------------------
__cleanup() {
    echo "Cleaning up intermediate directories..."
    rm -rf "${WORK_DIR}"/{process,calibrated,light_c,light_c_a}
    echo "Cleanup complete."
}

__usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  --hot-pixel-map          Create and use bad pixel map to detect hot and cold pixels (and not master dark)
  --align                  Perform plate solving and alignment
  --mbias <file>           Use existing master bias file
  --mdark <file>           Use existing master dark file
  --mflat <file>           Use existing master flat file
  --cleanup                Remove intermediate processing directories
  --help                   Show this help message

Examples:
  $(basename "$0") --align
  $(basename "$0") --mbias mbias.fit --mdark mdark.fit --mflat mflat.fit
  $(basename "$0") --cleanup
EOF
    exit 0
}

#--------------------------------------------
# Variables and argument parsing.
#--------------------------------------------
OPTION_ALIGN=false
OPTION_HOTPIXEL_MAP=false
MBIAS_PATH=""
MDARK_PATH=""
MFLAT_PATH=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --cleanup)
      __cleanup
      exit 0
      ;;
    --help)
      __usage
      ;;
    --align)
      OPTION_ALIGN=true
      shift
      ;;
    --hot-pixel-map)
      OPTION_HOTPIXEL_MAP=true
      shift
      ;;
    --mbias)
      MBIAS_PATH="$2"
      shift 2
      ;;
    --mdark)
      MDARK_PATH="$2"
      shift 2
      ;;
    --mflat)
      MFLAT_PATH="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information."
      exit 1
      ;;
  esac
done

#--------------------------------------------
# Siril cli executable and work directory path.
#--------------------------------------------
WORK_DIR=$(pwd)
SIRIL_PATH=/Applications/Siril.app/Contents/MacOS/siril-cli

__cleanup

#--------------------------------------------
# Setup Siril script for additional concat.
#--------------------------------------------
SIRIL_SCRIPT="requires 1.0.0
setext fits
set16bits
"

#--------------------------------------------
# Bias.
#--------------------------------------------
if [ -z "$MBIAS_PATH" ]; then
    SIRIL_SCRIPT+="
cd bias
convert bias -out=../process
cd ../process
stack bias rej 3 3 -nonorm -out=../master/mbias
cd ..
"
else
    cp "$MBIAS_PATH" master/mbias.fits
    SIRIL_SCRIPT+="# Using provided master bias: $MBIAS_PATH
"
fi

#--------------------------------------------
# Flat.
#--------------------------------------------
if [ -z "$MFLAT_PATH" ]; then
    SIRIL_SCRIPT+="
cd flat
convert flat -out=../process
cd ../process

# Calibrate flat frames by subtracting mbias
calibrate flat -bias=../master/mbias

# Stack flat frames to mflat.fit
stack pp_flat rej 3 3 -norm=mul -out=../master/mflat
cd ..
"
else
    cp "$MFLAT_PATH" master/mflat.fits
    SIRIL_SCRIPT+="# Using provided master flat: $MFLAT_PATH
"
fi

#--------------------------------------------
# Dark.
#--------------------------------------------
if [ -z "$MDARK_PATH" ]; then
    SIRIL_SCRIPT+="
cd dark
convert dark -out=../process
cd ../process

# Stack dark frames to mdark.fit
stack dark rej 3 3 -nonorm -out=../master/mdark
cd ..
"
else
    cp "$MDARK_PATH" master/mdark.fits
    SIRIL_SCRIPT+="# Using provided master dark: $MDARK_PATH
"
fi

#--------------------------------------------
# Create bad pixel map by creating median dark.
#--------------------------------------------
if [ "$OPTION_HOTPIXEL_MAP" = true ]; then
    if [ ! -f "./pixels.lst" ]; then
	SIRIL_SCRIPT+="
cd dark
convert dark -out=../process
cd ../process

# Median master stack for hot/cold pixel finding.
stack dark median -nonorm -out=../master/mmediandark
cd ..

# Load median master dark and create find hot/cold pixel list file.
load master/mmediandark
find_hot pixels 3 3
"
    else
	SIRIL_SCRIPT+="
# Using existing bad pixel map file pixels.lst
"
    fi
fi


#--------------------------------------------
# Light.
#--------------------------------------------
SIRIL_SCRIPT+="
cd light
convert light -out=../calibrated
cd ../calibrated
"

#--------------------------------------------
# Hot pixel map.
#--------------------------------------------
# Calibrate and detect hot and cold pixels from the masterdark.
if [ "$OPTION_HOTPIXEL_MAP" = false ]; then
    SIRIL_SCRIPT+="
calibrate light -dark=../master/mdark -flat=../master/mflat -cc=dark
"
else # Calibrate and with bad pixel map created with find_hot.
    SIRIL_SCRIPT+="
calibrate light -dark=../master/mdark -flat=../master/mflat -cc=bpm ../pixels.lst
"
fi

#--------------------------------------------
# WCS alignment.
#--------------------------------------------
if [ "$OPTION_ALIGN" = true ]; then
    SIRIL_SCRIPT+="

load pp_light_00001
parse '$RA:ra$_$DEC:dec$'
platesolve -force -disto=platesolve_data.wcs

register pp_light -disto=file platesolve_data.wcs
"    
fi

SIRIL_SCRIPT+="close"

echo "$SIRIL_SCRIPT" | $SIRIL_PATH -d "$WORK_DIR" -s -

# Copy calibrated and WCS aligned files in to Tycho naming convention directories.
mkdir -p light_c && cp calibrated/pp_light*.fits light_c
if [ "$OPTION_ALIGN" = true ]; then
    mkdir -p light_c_a && cp calibrated/r_pp_light*.fits light_c_a
fi
