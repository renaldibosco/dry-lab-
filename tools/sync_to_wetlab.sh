#!/bin/sh
# Dry Lab is the source of truth for the files both apps share (shared/assets/www/shared: base.css, i18n.js, ui.js, fonts).
# After changing them, copy them into the wet-lab checkout:   sh tools/sync_to_wetlab.sh ../wet-lab
set -e
DEST="${1:-../wet-lab}"
[ -d "$DEST/shared" ] || { echo "Cannot find $DEST/shared"; exit 1; }
rm -rf "$DEST/shared/assets/www/shared"
mkdir -p "$DEST/shared/assets/www"
cp -r shared/assets/www/shared "$DEST/shared/assets/www/shared"
echo "Copied shared files to $DEST"
