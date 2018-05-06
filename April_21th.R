?contour
#1st traintest result's follow-up analysis of its ols-cross validation matrix
cv.index.matrix = cv.list$ols.cv
contour(x = exp(seq(log(0.1),log(40), length=40)), y =seq(0,40, length=41), z = cv.index.matrix,
        xlim = c(0.1, 40), ylim = c(0,40),zlim = range(cv.index.matrix), levels = pretty(range(cv.index.matrix), 1000))

which.max(cv2$Freq)
cv2[1330:1350,]
range(cv.index.matrix)

#the ols is very unstable, so the smoothing when plot the contour will ignore some points while only display the obvious pattern