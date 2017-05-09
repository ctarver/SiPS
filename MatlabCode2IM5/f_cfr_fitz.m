function v = f_cfr_fitz(z,L,N,Fs)

R=zeros(1,N);

for m=1:N,
   summa=0;
   for k=m+1:L,
      summa=summa+z(k)*conj(z(k-m));
   end
   R(m)=1/(L-m)*summa;
end

v=1/(pi*N*(N+1)/Fs)*sum(atan2(imag(R),real(R)));
