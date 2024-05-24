cp constants.hil python/hexedpy/constants.py
cd python
sed -i "s/\\^/**/g" hexedpy/constants.py
python3 -m build
