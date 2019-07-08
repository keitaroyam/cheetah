setup_script=/home/sacla_sfx_app/setup.sh
python=/home/sacla_sfx_app/packages/rayonix/dials-v1-10-2/build/bin/dials.python
cheetah=/home/sacla_sfx_app/packages/rayonix/cheetah-sacla-api-rayonix
crystfel=/home/sacla_sfx_app/local/bin
scr=/home/sacla_sfx_app/packages/tools

sed -e "s,@@SETUP_SCRIPT@@,$setup_script,; s,@@PYTHON@@,$python,; s,@@CHEETAH_PATH@@,$cheetah,; s,@@INDEXAMAJIG_PATH@@,$crystfel,; s,@@SCRIPT_PATH@@,$scr," cheetah_dispatcher.py > cheetah_dispatcher_jul2019.py
