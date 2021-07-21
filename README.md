# circadian_shiny
## change log

### 2021/03/27
增添了批处理绘制热图，请前往**batch tab**页。重构了`rep_reduce`函数提取子dataframe。添加了默认图片的缓存。

- [x] 采用**Shiny Modules**拆分重构

- [x] 添加不同基因组id转换功能

- [ ] 尝试采用[*redis*](https://redis.io/)进行disk cache, [Rstudio ref article](https://shiny.rstudio.com/articles/caching.html)

### 2021/04/05
模块化拆分基本完成，添加了单细胞数据图像。

- [x] 优化单细胞测序结果模块，细化绘图逻辑

### 2021/04/09
添加了单细胞PCA降维图像，为单细胞模块建立了css文件，美化布局

- [ ] 为批处理提供更多功能

### 2021/04/12
设置了服务器每日重启，为不同组织的表达量水平提供了热图用于批处理

- [ ] 修复批处理需要resize界面的奇怪bug，大概率是与响应式相关，暂时想不到解决方案。

### 2021/7/21
新增了tool页面用于提供不同版本基因组之间的ID转换，暂时只提供G. hirsutum 中**WHU v1**基因组转**HAU v1**以及**ZJU v2.1**转**HAU v1**的功能。
- [ ] 提供更多基因组版本以供选择
- [ ] 整合不同版本基因组与拟南芥同于比对的数据