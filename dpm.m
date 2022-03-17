% 利用分布式幂方法估计Us矩阵。
% 初始化Us
Us = [];
% 设置独立信源的个数，文献1-6称为K
L = 2;
% 初始化各子阵的阵元个数，文献1-6称为Mk
antenna_in_subarray = 5;
% 设置子阵的个数
number_of_subarray = 3;
% 计算求得总阵元个数M，用于l的迭代
M = number_of_subarray * antenna_in_subarray;
% 设置一个足够大的AC迭代次数Q
Q = 10;