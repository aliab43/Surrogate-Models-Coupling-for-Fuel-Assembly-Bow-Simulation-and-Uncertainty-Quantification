import os
import sys
import subprocess
from glob import glob

def create_video_ffmpeg(frames_dir, out_file, total_seconds=60):
    frames = sorted(glob(os.path.join(frames_dir, 'frame_*.png')))
    if not frames:
        raise FileNotFoundError('No frames found in ' + frames_dir)
    n = len(frames)
    fps = max(1, min(60, int(round(n / total_seconds))))
    # ffmpeg expects an input pattern; ensure frame numbering is contiguous
    # Use image2 input with frame_%04d.png pattern
    pattern = os.path.join(frames_dir, 'frame_%04d.png')
    cmd = [
        'ffmpeg', '-y', '-framerate', str(fps), '-i', pattern,
        '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', out_file
    ]
    print('Running ffmpeg:', ' '.join(cmd))
    subprocess.check_call(cmd)
    print('Video written to', out_file)


def create_video_moviepy(frames_dir, out_file, total_seconds=60):
    from moviepy.editor import ImageSequenceClip
    frames = sorted(glob(os.path.join(frames_dir, 'frame_*.png')))
    if not frames:
        raise FileNotFoundError('No frames found in ' + frames_dir)
    n = len(frames)
    fps = max(1, min(60, float(n) / total_seconds))
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(out_file, codec='libx264')


def main():
    if len(sys.argv) < 2:
        print('Usage: python create_video.py <frames_dir> [out_file] [max_seconds]')
        sys.exit(1)
    frames_dir = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else 'deformation_video.mp4'
    total_seconds = int(sys.argv[3]) if len(sys.argv) > 3 else 60

    # Try ffmpeg first
    try:
        create_video_ffmpeg(frames_dir, out_file, total_seconds=total_seconds)
    except Exception as e:
        print('ffmpeg failed or not found, falling back to moviepy:', e)
        try:
            create_video_moviepy(frames_dir, out_file, total_seconds=total_seconds)
        except Exception as e2:
            print('moviepy fallback failed:', e2)
            raise

if __name__ == '__main__':
    main()
