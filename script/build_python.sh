cp constants.hil python/hexedpy/constants.py
cd python
sed -i "s/\\^/**/g" hexedpy/constants.py
sed -i "1i## \\\\namespace hexedpy.constants" hexedpy/constants.py
sed -i "1a# \\\\brief Simply ports \`hexed::constants\` into Python" hexedpy/constants.py
sed -i "2a# \\\\see \`hexed::constants\` for additional information" hexedpy/constants.py
python3 -m build
