function plot_sequence(seq)

len = size(seq, 2);
pos = -1:2/(len-1):1;
A_subseq = (seq=='A');
G_subseq = (seq=='G');
C_subseq = (seq=='C');
T_subseq = (seq=='T');

subplot(411);  stem(pos, A_subseq); axis([-1 1 0 2.0]); ylabel('A');
set(gca,'FontSize',18);
hold on
subplot(412);  stem(pos, G_subseq); axis([-1 1 0 2.0]); ylabel('G');
set(gca,'FontSize',18);
hold on 
subplot(413);  stem(pos, C_subseq); axis([-1 1 0 2.0]); ylabel('C');
set(gca,'FontSize',18);
hold on 
subplot(414);  stem(pos, T_subseq); axis([-1 1 0 2.0]); ylabel('T');
set(gca,'FontSize',18);
hold on 

end