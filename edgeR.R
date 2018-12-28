edgeR = function(data, param_G1, param_G2, filename){
  out_f1 <- sprintf("%s_edgeR.txt",filename)                  #出力ファイル名を指定してout_f1に格納
  out_f2 <- sprintf("%s_edgeR.png",filename)                  #出力ファイル名を指定してout_f2に格納
  param_FDR <- 0.05                      #false discovery rate (FDR)閾値を指定
  param_fig <- c(900, 900)               #ファイル出力時の横幅と縦幅を指定(単位はピクセル)
#必要なパッケージをロード
library(edgeR)                         #パッケージの読み込み
data.cl <- c(rep(1, param_G1), rep(2, param_G2))#G1群を1、G2群を2としたベクトルdata.clを作成

#edgeRを用いてDEGの検出を実行
d <- DGEList(counts=data,group=data.cl)#DGEListオブジェクトを作成してdに格納
d <- calcNormFactors(d)                #TMM正規化を実行
d <- estimateCommonDisp(d)             #the quantile-adjusted conditional maximum likelihood (qCML)法でcommon dispersionを計算している
d <- estimateTagwiseDisp(d)            #the quantile-adjusted conditional maximum likelihood (qCML)法でmoderated tagwise dispersionを計算している
out <- exactTest(d)                    #exact test (正確確率検定)で発現変動遺伝子を計算した結果をoutに格納
#tmp <- topTags(out, n=nrow(data), sort.by="none")#検定結果を抽出
p.value <- out$table$PValue            #p値をp.valueに格納
q.value <- p.adjust(p.value, method="BH")#q値をq.valueに格納
ranking <- rank(p.value)               #p.valueでランキングした結果をrankingに格納
# return(sum(q.value < param_FDR))               #FDR閾値(q.value < param_FDR)を満たす遺伝子数を表示

#ファイルに保存(テキストファイル)
tmp = as.data.frame(cbind(data, p.value, q.value, ranking)) #入力データの右側にp.value、q.value、rankingを結合した結果をtmpに格納
rownames(tmp) = rownames(data)
colnames(tmp) = colnames(tmp)
return(tmp)
# ed  
# #RPM正規化後のデータでM-A plotを描画（するための基礎情報取得）
# d <- DGEList(counts=data,group=data.cl)#DGEListオブジェクトを作成してdに格納
# d <- calcNormFactors(d)                #TMM正規化係数を計算
# norm_f_TMM <- d$samples$norm.factors   #TMM正規化係数の情報を抽出してnorm_f_TMMに格納
# names(norm_f_TMM) <- colnames(data)    #norm_f_TMMのnames属性をcolnames(data)で与えている
# effective_libsizes <- colSums(data) * norm_f_TMM#effective library sizesというのはlibrary sizesに(TMM)正規化係数を掛けたものなのでそれを計算した結果をeffective_libsizesに格納
# RPM_TMM <- sweep(data, 2, 1000000/effective_libsizes, "*")#元のカウントデータをeffective_libsizesで割り（RPMデータと同程度の数値分布にしたいので）1000000を掛けた正規化後のデータをRPM_TMMに格納
# 
# data <- RPM_TMM                        #RPM_TMMをdataに格納
# mean_G1 <- log2(apply(as.matrix(data[,data.cl==1]), 1, mean))#遺伝子ごとにG1群の平均の対数を計算した結果をmean_G1に格納
# mean_G2 <- log2(apply(as.matrix(data[,data.cl==2]), 1, mean))#遺伝子ごとにG2群の平均の対数を計算した結果をmean_G2に格納
# x_axis <- (mean_G1 + mean_G2)/2        #「G1群の平均値」と「G2群の平均値」の平均をとったものがM-A plotのA(x軸の値)に相当するものなのでx_axisに格納)
# y_axis <- mean_G2 - mean_G1            #いわゆるlog比(logの世界での引き算)がM-A plotのM(y軸の値)に相当するものなのでy_axisに格納)
# DEG_posi <- (stat_edgeR < param_FDR)      #指定した閾値未満のものの位置情報をDEG_posiに格納
# 
# #MA-plotを描画（本番）
# png(out_f2, pointsize=13, width=param_fig[1], height=param_fig[2])#出力ファイルの各種パラメータを指定
# plot(x_axis, y_axis, xlab="A=(log2(G2)+log2(G1))/2", ylab="M=log2(G2)-log2(G1)", pch=20, cex=.1)#MA-plotを描画
# grid(col="gray", lty="dotted")         #指定したパラメータでグリッドを表示
# points(x_axis[DEG_posi], y_axis[DEG_posi], col="magenta", pch=20, cex=0.1)#DEGを赤色にしている
# dev.off()                              #おまじない
# sum(DEG_posi)                          #発現変動遺伝子数を表示
# sum(DEG_posi)/nrow(data)               #発現変動遺伝子の全遺伝子数に占める割合を表示
}