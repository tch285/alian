#!/bin/bash

function abspath()
{
  case "${1}" in
    [./]*)
    echo "$(cd ${1%/*}; pwd)/${1##*/}"
    ;;
    *)
    echo "${PWD}/${1}"
    ;;
  esac
}
export -f abspath

# echo "This is only doing symlinks to yasp.."

echo "[i] Recipe dir is: {{yasp.recipe_dir}}"
#yasp --exec alian_dir=dirname {{yasp.recipe_dir}}

module_name=alian

version=current

rm -rf {{prefix}}/lib/alian
mkdir -p {{prefix}}/lib
ln -sv {{alian_dir}}/alian {{prefix}}/lib

mkdir -p {{prefix}}/bin

recipe_file={{yasp.recipe_file}}
this_recipe_name={{yasp.recipe}}
current_dir={{yasp.current_dir}}

shower_flag="-DCUSTOM_SHOWER=OFF"
if [ "{{custom_shower}}" = "on" ]; then
    shower_flag=-DCUSTOM_SHOWER=ON
    echo "[i] custom shower option is: ${shower_flag} -- called with custom_shower={{custom_shower}}"
else
  echo "[i] custom shower option is: ${shower_flag} -- use '--opt custom_shower=on' to enable it -- called with custom_shower={{custom_shower}}"
fi

build_visual="-DBUILD_VISUAL=OFF"
if [ "{{visual}}" = "on" ]; then
    build_visual=-DBUILD_VISUAL=ON
    echo "[i] custom shower option is: ${build_visual} -- called with visual={{visual}}"
else
  echo "[i] custom shower option is: ${build_visual} -- use '--opt visual=on' to enable it -- called with visual={{visual}}"
fi

srcdir=$(abspath {{yasp.recipe_dir}}/../alian)
cd {{builddir}}
echo "[i] source dir is {{srcdir}}"
cmake -DCMAKE_INSTALL_PREFIX={{prefix}} \
	-DCMAKE_BUILD_TYPE=Release ${shower_flag} ${build_visual} \
	{{srcdir}} && cmake --build . --target install -- -j {{n_cores}}
# {{alian_dir}}
