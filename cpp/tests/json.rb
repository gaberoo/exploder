
require 'json'
require 'pp'

n = 100
lamb = Array.new(n+1)
mu = Array.new(n+1)
psi = Array.new(n+1)

(0..n).to_a.each do |i|
  lamb[i] = 1.0 - i/n.to_f
  mu[i] = 0.1
  psi[i] = 0.1
end

out = Hash.new
out["lambda"] = lamb
out["mu"] = mu
out["psi"] = psi

puts out.to_json

