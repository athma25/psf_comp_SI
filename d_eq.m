function out=d_eq(n,pop,s,sp,sg_s,alp,et,nu,et0,s0,imm)

	pop_d=pop.*(1-((sp-s).^4./(sg_s.^4))-nu.*alp*pop)+imm;
	s_d=et0*(s0-s)+(et.*(sp-s))'*pop;
	out=[pop_d; s_d];
end
