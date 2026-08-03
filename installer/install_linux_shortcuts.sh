#!/bin/sh

set -eu

appimage=${APPIMAGE:-}
appdir=${APPDIR:-}
home=${HOME:-}

if [ -z "$appimage" ] || [ -z "$appdir" ] || [ -z "$home" ] || [ ! -f "$appimage" ]; then
    exit 0
fi

data_home=${XDG_DATA_HOME:-"$home/.local/share"}
applications_dir="$data_home/applications"
icon_dir="$data_home/icons/hicolor/256x256/apps"
launcher="$applications_dir/hip2linteractomics.desktop"
icon_source="$appdir/usr/share/icons/hicolor/256x256/apps/hip2linteractomics.png"

mkdir -p "$applications_dir" "$icon_dir"
if [ -f "$icon_source" ]; then
    cp -f "$icon_source" "$icon_dir/hip2linteractomics.png"
    chmod 644 "$icon_dir/hip2linteractomics.png"
fi

escaped_appimage=$(printf '%s' "$appimage" | sed 's/\\/\\\\/g; s/"/\\"/g; s/`/\\`/g; s/\$/\\$/g')
temporary="$launcher.tmp.$$"
umask 022
{
    printf '%s\n' '[Desktop Entry]'
    printf '%s\n' 'Version=1.0'
    printf '%s\n' 'Type=Application'
    printf '%s\n' 'Name=HIP2LInterActomics'
    printf '%s\n' 'Comment=Intermolecular interaction analysis with LUNA'
    printf 'Exec="%s"\n' "$escaped_appimage"
    printf '%s\n' 'Icon=hip2linteractomics'
    printf '%s\n' 'Categories=Science;Chemistry;'
    printf '%s\n' 'Terminal=false'
    printf '%s\n' 'StartupNotify=true'
} > "$temporary"
mv -f "$temporary" "$launcher"
chmod 644 "$launcher"

desktop_dir=""
if command -v xdg-user-dir >/dev/null 2>&1; then
    desktop_dir=$(xdg-user-dir DESKTOP 2>/dev/null || true)
fi
if [ -z "$desktop_dir" ] || [ "$desktop_dir" = "$home" ]; then
    desktop_dir="$home/Desktop"
fi
mkdir -p "$desktop_dir"
desktop_shortcut="$desktop_dir/HIP2LInterActomics.desktop"
cp -f "$launcher" "$desktop_shortcut"
chmod 755 "$desktop_shortcut"

if command -v gio >/dev/null 2>&1; then
    gio set "$desktop_shortcut" metadata::trusted true >/dev/null 2>&1 || true
fi
if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database "$applications_dir" >/dev/null 2>&1 || true
fi

bin_dir="$home/.local/bin"
mkdir -p "$bin_dir"
for command_name in hipplinteractomics-terminal hipplinteractomics-multiple-run; do
    command_flag="--terminal"
    if [ "$command_name" = "hipplinteractomics-multiple-run" ]; then
        command_flag="--multiple-run"
    fi
    command_path="$bin_dir/$command_name"
    command_tmp="$command_path.tmp.$$"
    {
        printf '%s\n' '#!/bin/sh'
        printf 'exec "%s" "%s" "$@"\n' "$escaped_appimage" "$command_flag"
    } > "$command_tmp"
    mv -f "$command_tmp" "$command_path"
    chmod 755 "$command_path"
done

profile="$home/.profile"
path_line='export PATH="$HOME/.local/bin:$PATH"'
if ! grep -F "$path_line" "$profile" >/dev/null 2>&1; then
    printf '\n%s\n' "$path_line" >> "$profile"
fi
