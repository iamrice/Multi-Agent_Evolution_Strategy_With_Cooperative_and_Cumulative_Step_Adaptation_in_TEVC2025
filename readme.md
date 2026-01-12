
## Usage

When test on benchmark functions：
```
cd ./main
g++ -std=c++11 CCSA_DES.cpp ../Benchmarks/Benchmarks.cpp ../Benchmarks/WSN_Benchmarks.cpp -o test
./test F1 (You can replace ‘F1‘ by other functions，just refer to the file “./Benchmarks/default_config.json“)
```

When test on application problems：

```
cd ./main
g++ -std=c++11 CCSA_DES_wsn.cpp ../Benchmarks/WSN_Benchmarks.cpp ../Benchmarks/Benchmarks.cpp -o wsn_task
./wsn_task WSN_location_RSS_3d_16s_10t (You can replace it by other application problems，just refer to the file ‘./Benchmarks/default_config.json‘)
```


## Reference

[T. -Y. Chen, W. -N. Chen, J. -K. Hao, Y. Wang and J. Zhang, "Multiagent Evolution Strategy With Cooperative and Cumulative Step Adaptation for Black-Box Distributed Optimization," in IEEE Transactions on Evolutionary Computation, vol. 29, no. 6, pp. 2819-2833, Dec. 2025, doi: 10.1109/TEVC.2025.3525713.](https://ieeexplore.ieee.org/document/10824905)

## Authors

**Wei-Neng Chen** (Senior Member, IEEE)
South China University of Technology, China.
Email: cschenwn@scut.edu.cn

**Wei-Neng Chen** received the bachelor’s and Ph.D. degrees in computer science from Sun Yat-sen University, Guangzhou, China, in 2006 and 2012, respectively. Since 2016, he has been a Full Professor with the School of Computer Science and Engineering, South China University of Technology, Guangzhou. He has coauthored over 100 international journal and conference papers, including more than 70 papers published in the IEEE TRANSACTIONS journals. His current research interests include computational intelligence, swarm intelligence, network science, and their applications. Dr. Chen was a recipient of the IEEE Computational Intelligence Society Outstanding Dissertation Award in 2016 and the National Science Fund for Excellent Young Scholars in 2016. He was also a PI of the National Science and Technology Innovation 2030—the Next Generation Artificial Intelligence Key Project. He is currently the Vice-Chair of the IEEE Guangzhou Section, and the Chair of IEEE SMC Society Guangzhou Chapter. He is also a Committee Member of the IEEE CIS Emerging Topics Task Force. He serves as an Associate Editor for the IEEE TRANSACTIONS ON NEURAL NETWORKS AND LEARNING SYSTEMS and the Complex and Intelligent Systems.

**Tai-You Chen** (Student Member, IEEE)
South China University of Technology, China.
Email: cstaiutan@mail.scut.edu.cn

**Tai-You Chen** received the bachelor’s degree in computer science and technology from the South China University of Technology, Guangzhou, China, in 2022, where he is currently pursuing the Ph.D. degree in computer science and technology with the School of Computer Science and Engineering.  His current research interests include swarm intelligence, evolutionary computation, consensus-based distributed optimization, multi-agent systems, and their applications in real-world problems.
