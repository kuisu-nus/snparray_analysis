# snparray_analysis
用于处理Genome Studio的SNP array数据，并将其转为VCF数据


### 匹配Experiment的异常情况记录
Cyto-12,R-温善荣,T11,T12,T13,T14,T15,T16,李笑玲B  L1,L2,L3        （未完）,刘美L4（重上）,RB,

```markdown
Cyto-12,F-潘毅,陈晓燕H   C1,C2,C3,詹杜鹃Z1,Z2,Z3,Z4,Z5,Z6,Z7,RB,
Cyto-12,M-陈黄花,潘燕婷P1,P2,P3,P4,李惠娴L1,L2,L3,L4,L5,L6,RB,
Cyto-12,R-潘丽婷,L7,L8,L9,L10,L11,L12,王生兰    R1,林木秀   R1,质控1,质控2,RB,
,R01C01,R02C01,R03C01,R04C01,R05C01,R06C01,R01C02,R02C02,R03C02,R04C02,R05C02,R06C02,
Cyto-12,M-楼英玲,曾佩环Z1,Z2,Z3,Z4,江丽红B   J1,J2,J3,J4,J5,J6,RB,
Karyomap,F-叶久思,M-周玲娟(预实验),R-叶檬,R-董依依,郑乐美Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,
Karyomap,Z9,Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,
,X1 number：203218670040                                ,,,,,,X2 number：203218670041,,,,,,
,X3 number：203359230010                               ,,,,,,X4 number：203359230011,,,,,,
,X5 number：203392440048                               ,,,,,,X6 number：203392440054,,,,,,
```
这里一个实验记录了6张芯片，需要修改策略

```markdown
,,,,,,,,,,,,,
,实验时间： 2019  年 10 月 29  日,,,,,,数据分析：  潘家富                                                    打印第   页,,,,,,
Cyto-12,F-何卓良,卢洁红L1,L5,"陈晓梅C
C3","林琳H
L1",L2,何雪娇H1,H2,H3,H4,H7,RB,
Cyto-12,M-司徒凤仪,蔡晓欣C1,司徒凤仪S1,S2,S3,S4,S5,S6,S7,S8,S9,RB,
,R01C01,R02C01,R03C01,R04C01,R05C01,R06C01,R01C02,R02C02,R03C02,R04C02,R05C02,R06C02,
Cyto-12,R-梁凤娥,陈倩霞C1,C2,C3,C4,C5,C6,C7,C8,C9,杨丽G1,RB,
Karyomap,F-李振亚,M-司秋锦(预实验),R-司秋锦儿子,邹承娇Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,
,X1 number：203880360095                                ,,,,,,X2 number：203880360099,,,,,,
,X3 number：203880360128                               ,,,,,,X4 number：203717980048,,,,,,
备注： 1.加粗字体为二次活检样本。,,,,,,,,,,,,,
备注：红色字体表示需要区分正常和携带；Cyto-12芯片中R01C01为阳性对照外周血，如无特殊备注则为同一家系外周血：F表示父亲，M-母亲，R-先证者,,,,,,,,,,,,,
           标绿背景的病人表示按新收费标准,,,,,,,,,,,,,
,,,,,,,,,,,,,
````
这里，每个一行中有多个换行符


```markdown
芯片实验记录表,,,,,,,,,,,,,
,实验时间：2020  年  月   日,,,,,,数据分析：                                                         打印第   页,,,,,,
Cyto-12,F-,,,,,,,,,,,RB,
Cyto-12,M- ,,,,,,,,,,,RB,
,R01C01,R02C01,R03C01,R04C01,R05C01,R06C01,R01C02,R02C02,R03C02,R04C02,R05C02,R06C02,
Cyto-12,R-,,,,,,,,,,,RB,
Karyomap,F-,M-(预实验),R-,,,,,,,,,,
,X1 number：                                ,,,,,,X2 number：,,,,,,
,X3 number：                               ,,,,,,X4 number：,,,,,,
备注：红色字体表示需要区分正常和携带；Cyto-12芯片中R01C01为阳性对照外周血，如无特殊备注则为同一家系外周血：F表示父亲，M-母亲，R-先证者,,,,,,,,,,,,,
           标绿背景的病人表示按新收费标准,,,,,,,,,,,,,
```
这个数据中没有任何信息

```markdown
025-10-27 16:10:40,662 - WARNING - Unable to parse name:
2025-10-27 16:10:40,662 - WARNING - Unable to parse name:
2025-10-27 16:10:40,662 - WARNING - Unable to parse name:
2025-10-27 16:10:40,663 - WARNING - Unable to parse name: 刘秋连B L6？
2025-10-27 16:10:40,663 - WARNING - Unable to parse name: 426
2025-10-27 16:10:40,663 - WARNING - Unable to parse name: 106NC
2025-10-27 16:10:40,664 - WARNING - Unable to parse name: 12t
2025-10-27 16:10:40,666 - WARNING - Unable to parse name: 聂红群（羊水）
2025-10-27 16:10:40,667 - WARNING - Unable to parse name: W21
2025-10-27 16:10:40,667 - WARNING - Unable to parse name: 4
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V2 1+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V2 4+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V2 10+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 容燕玲V2 1+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 蔡丽琳V2  1
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 钟丽娟V2  2
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V3 13+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V3 14+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 韦玉美V3 17+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 容燕玲V3  2+
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 蔡丽琳V3  3
2025-10-27 16:10:40,669 - WARNING - Unable to parse name: 蔡丽琳V3  4
2025-10-27 16:10:40,670 - WARNING - Unable to parse name: 廖永丽   (1、1)
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖永丽 (1、2)
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖永丽 (1、3)
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖永丽 (1、4)
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖永丽 (1、5)
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖永丽B
2025-10-27 16:10:40,671 - WARNING - Unable to parse name: 廖兰秀   (2、1)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 廖兰秀   (2、2)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 廖兰秀  (2、3)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 廖兰秀    (2、4)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 廖兰秀    (2、5)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 廖兰秀B
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 卢楚韵   (3、1)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 卢楚韵   (3、2)
2025-10-27 16:10:40,672 - WARNING - Unable to parse name: 卢楚韵   (3、3)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 卢楚韵   (3、4)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 卢楚韵   (3、5)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 卢楚韵B (脐血)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 张玲P   (4.1)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 张玲P   (4.2)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 张玲P   (4.3)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 张玲P   (4.4)
2025-10-27 16:10:40,673 - WARNING - Unable to parse name: 张玲P   (4.5)
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 张玲P    (脐血)
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 欧阳绮雯(羊水)
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 康为世纪 (MDA7+)
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 8+
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 9+
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 10+
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: Picoplex Gold    2+
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 3+
2025-10-27 16:10:40,674 - WARNING - Unable to parse name: 5+
2025-10-27 16:10:40,676 - WARNING - Unable to parse name: 梁媛   (脐血)
2025-10-27 16:10:40,677 - WARNING - Unable to parse name: 康为世纪1
2025-10-27 16:10:40,677 - WARNING - Unable to parse name: 康为世纪2
2025-10-27 16:10:40,677 - WARNING - Unable to parse name: 康为世纪3
2025-10-27 16:10:40,677 - WARNING - Unable to parse name: 康为世纪4
2025-10-27 16:10:40,678 - WARNING - Unable to parse name: 6.1     欧阳绮文
2025-10-27 16:10:40,678 - WARNING - Unable to parse name: 6.2           欧阳绮文
2025-10-27 16:10:40,678 - WARNING - Unable to parse name: 6.3           欧阳绮文
2025-10-27 16:10:40,678 - WARNING - Unable to parse name: 6.4           欧阳绮文
2025-10-27 16:10:40,678 - WARNING - Unable to parse name: 6.5           欧阳绮文
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 6              欧阳绮文
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 6              张浩峰
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 6    欧阳绮文B
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 7.1邱素芳
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 7.2邱素芳
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 7.3邱素芳
2025-10-27 16:10:40,679 - WARNING - Unable to parse name: 7.4邱素芳
2025-10-27 16:10:40,680 - WARNING - Unable to parse name: 7.5邱素芳
2025-10-27 16:10:40,680 - WARNING - Unable to parse name: 邱素芳B (脐血)
2025-10-27 16:10:40,680 - WARNING - Unable to parse name: 谢志梅ＢX1
2025-10-27 16:10:40,682 - WARNING - Unable to parse name: 李松 (P1D)
2025-10-27 16:10:40,682 - WARNING - Unable to parse name: 李松  (Q2D)
2025-10-27 16:10:40,682 - WARNING - Unable to parse name: 李松  (Q3D)
2025-10-27 16:10:40,682 - WARNING - Unable to parse name: 刘晓丹  胎盘1
2025-10-27 16:10:40,683 - WARNING - Unable to parse name: 刘晓丹  胎盘2
2025-10-27 16:10:40,683 - WARNING - Unable to parse name: 刘晓丹  胎盘3
2025-10-27 16:10:40,683 - WARNING - Unable to parse name: 刘晓丹  胎盘4
2025-10-27 16:10:40,683 - WARNING - Unable to parse name: 刘晓丹  胎盘5
2025-10-27 16:10:40,684 - WARNING - Unable to parse name: 9
2025-10-27 16:10:40,687 - WARNING - Unable to parse name: 康为世纪6
2025-10-27 16:10:40,687 - WARNING - Unable to parse name: 7
2025-10-27 16:10:40,687 - WARNING - Unable to parse name: 8
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 黄飘飘  胎盘1
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 黄飘飘  胎盘2
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 黄飘飘  胎盘3
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 黄飘飘  胎盘4
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 黄飘飘  胎盘5
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 吴媚媚  胎盘1
2025-10-27 16:10:40,688 - WARNING - Unable to parse name: 吴媚媚  胎盘2
2025-10-27 16:10:40,689 - WARNING - Unable to parse name: 吴媚媚  胎盘3
2025-10-27 16:10:40,689 - WARNING - Unable to parse name: 吴媚媚  胎盘4
2025-10-27 16:10:40,689 - WARNING - Unable to parse name: 吴媚媚  胎盘5
```

异常标签的数据


## 数据清洗

- 发现同一个病人，可能存在多个周期。
- 所有要多个标签进行约束: name-PGD编号-年份-遗传材料接收者