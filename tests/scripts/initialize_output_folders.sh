OUTT="$1"
if [ ! -d "test_outputs" ]; then
  mkdir "test_outputs"
fi
if [ -n "$OUTT" ] && [ -d "$OUTT" ]; then
  rm -rf "$OUTT"
fi
