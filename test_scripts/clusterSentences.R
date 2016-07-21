
# imperfect attempt to auto-cluster samples by description
ignore.words <- c("","MANC", "UoW", "UCSD", "Broad", "F", "M", "U")
sentence.vec <- colnames(data.subset)
all.words <- strsplit(sentence.vec, "\\.")
unique.words <- unique(unlist(all.words))
allowed.words <- setdiff(unique.words, ignore.words)
# what to do about partial matches (e.g. if a word is a single letter.)
allowed.words <- allowed.words[order(nchar(allowed.words), decreasing=T)]  # order by longest
grep(allowed.words[1], sentence.vec)
max.words <- length(allowed.words)  # how many words to use
score.matrix <- matrix(0, ncol=length(sentence.vec), nrow=max.words, dimnames=list(allowed.words[1:max.words], make.names(sentence.vec)))
for(i in 1:max.words) {
  hits <- grep(allowed.words[i], sentence.vec)
  score.matrix[i,hits] <- score.matrix[i,hits] + 1
}
#hclust(score.matrix)

plot(hclust(dist(t(score.matrix))))
