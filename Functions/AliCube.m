function [i1,i2,i3,i4,i5,i6]=AliCube(c,c1,c2,c3)


a = -pi : pi/2 : pi;
ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
x = [cos(a+ph); cos(a+ph)]/cos(ph);
y = [sin(a+ph); sin(a+ph)]/sin(ph);
z = [-ones(size(a)); ones(size(a))];
[M, ~, ImageAlpha] = imread('Alicube.png');
xe = [1 -1;1 -1]*c; ye = [1 1;-1 -1]*c; ze = [0.68 0.68;0.68 0.68]*c;
i1=surf(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S1');
 hold on
 xe = [1 -1;1 -1]*c; ye = [0.675 0.675;0.675 0.675]*c; ze = [1 1;-1 -1]*c;
i2=mesh(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S2');
xe = [1 -1;1 -1]*c; ye = [1 1;-1 -1]*c; ze = [-0.64 -0.64;-0.64 -0.64]*c;
i3=mesh(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S3');
xe = [1 -1;1 -1]*c; ye = [-0.64 -0.64;-0.64 -0.64]*c; ze = [1 1;-1 -1]*c;
i4=mesh(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S4');
xe = [0.52 0.52;0.52 0.52]*c; ye = [1.268 -1.2;1.268 -1.2]*c; ze = [1 1;-1 -1]*c;
i5=mesh(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S5');
xe = [-0.548 -0.548;-0.548 -0.548]*c; ye = [1.268 -1.19;1.268 -1.19]*c; ze = [1 1;-1 -1]*c;
i6=mesh(xe+c1,ye+c2,ze+c3,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none','tag','S6');
end