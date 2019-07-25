kvh1000 = read_kvh('/home/spiels/log/2019MAGEXPS/kvh/2019_07_18_14_08.KVH');
kvh100 = read_kvh('/home/spiels/log/2019MAGEXPS/kvh/2019_07_18_14_13.KVH');
kvh10 = read_kvh('/home/spiels/log/2019MAGEXPS/kvh/2019_07_18_14_26.KVH');

sensor_report(kvh10,'~/kvh10.pdf');
sensor_report(kvh100,'~/kvh100.pdf');
sensor_report(kvh1000,'~/kvh1000.pdf');
! pdftk kvh*.pdf cat output Kvh_noise_charac.pdf
! rm kvh10*.pdf