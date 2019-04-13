#!/bin/sh

root=/home/oys/cheetah/eiger-zmq/boost_python
export LD_LIBRARY_PATH=$root:$LD_LIBRARY_PATH
export PYTHONPATH=$root:$PYTHONPATH
echo $PYTHONPATH

PHENIX_TRUST_OTHER_ENV=1 yamtbx.python $root/cheetah_client.py \
 --bl=32xu --nproc=4 --logdir=log \
 --vent-hosts="127.0.0.1:9999" \
 --result-host="127.0.0.1:5558" \
 --pub-host=""  &

yamtbx.shika_backend nproc=0 bl=32xu > log/shika_backend.log 2>&1 &

wait
