function [ new_S_ON, new_SOFF ] = S_ON_SOFF_processor(S_ON, SOFF)
% S_ON_SOFF_processor
%   realigns the SOFF and S_ON for static
S_ON_index=1;
SOFF_index=1;
new_S_ON_index=1;
new_SOFF_index=1;
S_ON_delta_index = 0;
SOFF_delta_index = 0;
% the strategy is to always find the next SOFF that matches an S_ON
% if there should be one or more S_ON's in between an S_ON and a SOFF then
% that new S_ON has to replace the earlier S_ON which has no SOFF partner
made_new_S_ON_SOFF = 0;
while S_ON_index <=length(S_ON) &&  SOFF_index <=length(SOFF)
  new_S_ON(new_S_ON_index) = S_ON(S_ON_index);
  new_SOFF(new_SOFF_index) = SOFF(SOFF_index);  
  % check to see if S_ON is greater than SOFF
  while S_ON_index+1<= length(S_ON) && S_ON(S_ON_index)>SOFF(SOFF_index)
      % if so then look for a larger S_OFF that would match this S_ON
      SOFF_index = SOFF_index + 1;
      new_SOFF(new_SOFF_index) = SOFF(SOFF_index);
      made_new_S_ON_SOFF = 1;
  end
  new_S_ON_index = new_S_ON_index + 1;
  new_SOFF_index = new_SOFF_index + 1;
  S_ON_index = S_ON_index + 1;
  SOFF_index = SOFF_index + 1;
end

if made_new_S_ON_SOFF
    stim_durs = new_SOFF-new_S_ON;
    figure
    title(['made new S_ON SOFFs of size ' num2str(length(new_S_ON)) ', min(stim_durs) = ' ...
        num2str(min(stim_durs)) ', max(stim_durs) = ', num2str(max(stim_durs)) ...
        ', mean(stim_durs)=' num2str(mean(stim_durs))],'fontsize', 12, 'Interpreter','none')
    hold on
    hist(stim_durs, 400)
end
end
