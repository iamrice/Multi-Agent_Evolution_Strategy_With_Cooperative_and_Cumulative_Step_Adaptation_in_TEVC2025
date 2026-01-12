```
cd ./main
g++ -std=c++11 CCSA_DES.cpp ../Benchmarks/Benchmarks.cpp ../Benchmarks/WSN_Benchmarks.cpp -o test
./test F1 (可替换为其他测试函数，参考./Benchmarks/default_config.json)
```

```
cd ./main
g++ -std=c++11 CCSA_DES_wsn.cpp ../Benchmarks/WSN_Benchmarks.cpp ../Benchmarks/Benchmarks.cpp -o wsn_task
./wsn_task WSN_location_RSS_3d_16s_10t (可替换为其他测试函数，参考./Benchmarks/default_config.json)
```
