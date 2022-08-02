QMLdir='/home/robert/robert2/CODE_CORECT/QuantumModelLib'

echo $QMLdir
# exit
fpm build --link-flag "-fopenmp $QMLdir/libpot.a"

#fpm run --flag "-fopenmp /robert2/complo/QuantumModelLib/libpot.a
