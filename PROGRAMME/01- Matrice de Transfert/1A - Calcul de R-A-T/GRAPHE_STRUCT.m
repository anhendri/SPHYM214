A = load('OUTPUT.dat');
figure
plot(A(:,1),A(:,2:4),'LineWidth',1);
set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',1)
ylim([0,100])
%xlim([100,800])
ylabel('$\%$','Interpreter','LaTeX','FontSize',20)
xlabel('$\lambda (nm)$', 'Interpreter','LaTeX','FontSize',20)
legend('Réflectance','Absorbance')