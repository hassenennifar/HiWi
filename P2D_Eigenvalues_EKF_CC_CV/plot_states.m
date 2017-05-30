
subplot(5,2,5)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$c_l$ / mol m$^{-3}$', 'Interpreter', 'latex')
y = x(P.idx_cl:P.nx:P.nx*P.nj, :);
y_est = x_est( P.idx_cl:P.nx:P.nx*P.nj, :);

hold on
cla
%     plot(h, y(:,k), '-ro')
plot(com.p*1e6, com.d9(k,:), 'b')
plot(h, y_est(:,k), '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,6)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$i_n$ / A m$^{-2}$', 'Interpreter', 'latex')
y = x(P.idx_jn:P.nx:P.nx*P.nj, :)*P.F;
y_est = x_est(P.idx_jn:P.nx:P.nx*P.nj, :)*P.F;

hold on
cla
% 	plot(h, y(:,k), '-ro')
plot(com.p(1,1:54)*1e6, com.d10(k,1:54), 'b')
plot(h(1:P.bnd_sep_neg), y_est(1:P.bnd_sep_neg,k), '--ok', 'MarkerFaceColor', 'k')
if k<=2
leg =legend('Actual', 'Estimate');
set(leg, 'Interpreter', 'latex')
end
%     subplot(3,6,4)
%     y = x(P.idx_jn:P.nx:P.nx*P.nj, :);
%     y_est = x_est(P.idx_jn:P.nx:P.nx*P.nj, :);

hold on
%     cla
% 	plot(h, y(:,k), '-ro')
plot(com.p(1,70:end)*1e6, com.d10(k,70:end), 'b')
plot(h(P.bnd_pos_sep:P.nj), y_est(P.bnd_pos_sep:P.nj,k), '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,7)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$\Phi_l$ / V', 'Interpreter', 'latex')
y = x(P.idx_phil:P.nx:P.nx*P.nj, :);
y_est = x_est(P.idx_phil:P.nx:P.nx*P.nj, :)-x_est(P.idx_phis, k);

hold on
cla
%     plot(h, y(:,k), '-ro')
plot(com.p*1e6, com.d4(k,:), 'b')
plot(h, y_est(:,k), '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,9)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$i_l$ / A m$^{-2}$', 'Interpreter', 'latex')
xlabel('Cell length /(\mum)', 'Interpreter', 'Tex')
y = x(P.idx_il:P.nx:P.nx*P.nj, :);
y_est = x_est(P.idx_il:P.nx:P.nx*P.nj, :);

hold on
cla
%     plot(h, y(:,k), '-ro')
plot(com.p*1e6, com.d13(k,:), 'b')
plot(h, y_est(:,k), '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,8)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$c_{ss}/c_{s,\max}$', 'Interpreter', 'latex')
y = x(P.idx_css:P.nx:P.nx*P.nj, :);
y_est = x_est(P.idx_css:P.nx:P.nx*P.nj, :);

hold on
cla
% 	plot(h(1:P.bnd_sep_neg), y(1:P.bnd_sep_neg,k)/P.csmax_neg, '-ro')
plot(com.p(1,1:54)*1e6, com.d7(k,1:54)/P.csmax_neg, 'b')
plot(h(1:P.bnd_sep_neg), y_est(1:P.bnd_sep_neg,k)/P.csmax_neg, '--ok', 'MarkerFaceColor', 'k')

%     subplot(3,6,10)
y = x(P.idx_css:P.nx:P.nx*P.nj, :);
y_est = x_est(P.idx_css:P.nx:P.nx*P.nj, :);

hold on
%     cla
%     plot(h(P.bnd_pos_sep:P.nj), y(P.bnd_pos_sep:P.nj,k)/P.csmax_pos, '-ro')
plot(com.p(70:end)*1e6, com.d7(k,70:end)/P.csmax_pos, 'b')
plot(h(P.bnd_pos_sep:P.nj), y_est(P.bnd_pos_sep:P.nj,k)/P.csmax_pos, '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,10)

set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('$\eta$ / mV', 'Interpreter', 'latex')
xlabel('Cell length / \mum', 'Interpreter', 'Tex')
y = overpot(x,P);
y_est = overpot(x_est,P_est);
hold on
cla
plot(com.p(1,1:54)*1e6, com.d5(k,1:54)*1e3, 'b')
plot(h(1:P.bnd_sep_neg), y_est(1:P.bnd_sep_neg,k), '--ok', 'MarkerFaceColor', 'k')

plot(com.p(70:end)*1e6, com.d5(k,70:end)*1e3, 'b')
plot(h(P.bnd_pos_sep:P.nj), y_est(P.bnd_pos_sep:P.nj,k), '--ok', 'MarkerFaceColor', 'k')

%%
subplot(5,2,1:4)
set(gca, 'FontName', 'CMU Serif', 'FontSize', 11)
ylabel('Terminal voltage / V', 'Interpreter', 'latex')
cla
plot(ts(2:k), ym(2:k), 'b', 'LineWidth', 1)
hold on
plot(ts(2:k), ym1(2:k), '.', 'color', [.6,.6,.6])
plot(ts(1:k), v_est(1:k), 'k--')
if k <=2
leg=legend('Actual', 'Measured', 'Estimated');
set(leg, 'Interpreter', 'latex')
end


pause(.01)
% 
% if mod(k,2)==0
% frame = getframe(1);
% writeVideo(vidObj, frame);
% im = frame2im(frame);
% 
% [imind,cm] = rgb2ind(im,256);
% 
% if k == 2
%     imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',10);  
% elseif k ==50
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',10);
% elseif k ==length(ts)-1
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',10);
% else
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.01);
% end
% end


