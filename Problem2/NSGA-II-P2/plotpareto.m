function   plotpareto(F,C)


C1=C(F{1},:);


plot(C1(:,1),C1(:,2),'ko','markerfacecolor','k','MarkerSize',6)
xlabel(' F 1 ');
ylabel(' F 2 ');



end