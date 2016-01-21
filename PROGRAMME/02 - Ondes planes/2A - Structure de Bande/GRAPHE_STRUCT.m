A = load('OUTPUT.dat');
figure
plot(A(:,1),A(:,2:11),'LineWidth',1);
nk = size(A,1)-1;
set(gca,'XTickLabel',{'X','U','L','\Gamma','X','W','K'},...
    'XTick',[0 nk/6 nk/3 nk/2 2*nk/3 5*nk/6 nk], ...
    'FontSize',20,'FontName','Times New Roman','LineWidth',1)
ylim([0,2])
ylabel('$\omega \frac{a}{2\pi c}$','Interpreter','LaTeX','FontSize',20)