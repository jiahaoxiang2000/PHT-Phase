function [d]=solve_d(dstiff, drhs, d, dfree, dfixed)
d(dfree) = dstiff(dfree,dfree)\(drhs(dfree)-dstiff(dfree,dfixed)*d(dfixed));
